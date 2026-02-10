from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import os
import threading
import time

import gemmi  # pip install gemmi

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"         # your *.cif coordinate files
DSSP_DIR = BASE / "dssp_out"                 # where we store mkdssp outputs
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")
LIBCIFPP_DATA_DIR = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\share\libcifpp")

WORKERS = 8

# Increase if needed for huge structures
MKDSSP_TIMEOUT_SEC = 900  # 15 minutes per file

BETA_COUNTS_TSV = OUT_DIR / "beta_counts.tsv"
FAIL_LOG = OUT_DIR / "dssp_failures.log"
ZEROCOUNT_LOG = OUT_DIR / "dssp_zero_counts.log"

# DSSP beta codes
BETA_SS = {"E", "B"}
# ====================================================


def run_cmd(cmd: list[str], env: dict, timeout_sec: int) -> tuple[int, str, str, bool]:
    """
    Returns (returncode, stdout, stderr, timed_out).
    On timeout: kills the process and returns timed_out=True.
    """
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    try:
        so, se = p.communicate(timeout=timeout_sec)
        return p.returncode, so, se, False
    except subprocess.TimeoutExpired:
        try:
            p.kill()
        except Exception:
            pass
        so, se = p.communicate()
        return 124, so or "", (se or "") + f"\n[TIMEOUT after {timeout_sec}s]", True


def run_mkdssp_adaptive(cif_path: Path, out_text_path: Path, out_cif_path: Path) -> tuple[bool, Path, str, str]:
    """
    Try classic DSSP text first.
    If mkdssp complains about old DSSP format, rerun to produce mmCIF output.
    Returns: (ok, output_path_used, mode, err_msg)
      mode in {"dssp", "mmcif"}
    """
    env = os.environ.copy()
    env["LIBCIFPP_DATA_DIR"] = str(LIBCIFPP_DATA_DIR)

    # 1) Try classic DSSP text output: mkdssp input output
    rc, so, se, to = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_text_path)], env, MKDSSP_TIMEOUT_SEC)
    if rc == 0 and out_text_path.exists():
        return True, out_text_path, "dssp", ""

    err = (se or "") + (so or "")
    # If it fails due to "old DSSP format", switch to mmCIF output
    if "old DSSP format" in err or "mmCIF format instead" in err:
        # 2a) Try explicit output-format if supported
        rc2, so2, se2, to2 = run_cmd(
            [str(MKDSSP_EXE), "--output-format=mmcif", str(cif_path), str(out_cif_path)],
            env,
            MKDSSP_TIMEOUT_SEC,
        )
        if rc2 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif", ""

        # 2b) Fallback: some builds infer from extension
        rc3, so3, se3, to3 = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_cif_path)], env, MKDSSP_TIMEOUT_SEC)
        if rc3 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif", ""

        return False, out_cif_path, "mmcif", (se2 or "") + (so2 or "") + (se3 or "") + (so3 or "")

    # Other failure (including timeout)
    return False, out_text_path, "dssp", err


def parse_dssp_text_beta_counts(dssp_path: Path) -> Dict[str, int]:
    """
    Classic DSSP text parsing (fixed width).
    chain at col 12 (index 11), SS at col 17 (index 16).
    """
    beta: Dict[str, int] = {}
    lines = dssp_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    start = None
    for idx, line in enumerate(lines):
        if line.startswith("  #  RESIDUE"):
            start = idx + 1
            break
    if start is None:
        return beta

    for line in lines[start:]:
        if len(line) < 17:
            continue
        chain = line[11].strip()
        ss = line[16].strip()
        if chain and ss in BETA_SS:
            beta[chain] = beta.get(chain, 0) + 1
    return beta


def _find_tag_index(tags_lower: list[str], want_exact: list[str], want_suffix: list[str]) -> int | None:
    """
    Try exact matches first; then suffix matches.
    tags_lower are already lowercased.
    """
    for w in want_exact:
        w = w.lower()
        if w in tags_lower:
            return tags_lower.index(w)
    for suf in want_suffix:
        suf = suf.lower()
        for i, t in enumerate(tags_lower):
            if t.endswith(suf):
                return i
    return None


def parse_mkdssp_mmcif_beta_counts(mmcif_path: Path) -> Dict[str, int]:
    """
    Robust mmCIF parser for mkdssp output.

    Primary strategy:
      - Look for loop tags starting with "_dssp_struct_summary."
      - Find chain id column (auth_asym_id/label_asym_id/asym_id)
      - Find secondary structure column (secondary_structure/ss/structure)
      - Count E/B per chain

    Fallback strategy:
      - If that category isn't found, do a broader search.
    """
    doc = gemmi.cif.read_file(str(mmcif_path))
    block = doc.sole_block()

    # ---------- Primary: official DSSP category ----------
    for loop in block.loops:
        tags = [t.lower() for t in loop.tags]
        if not any(t.startswith("_dssp_struct_summary.") for t in tags):
            continue

        chain_idx = _find_tag_index(
            tags,
            want_exact=[
                "_dssp_struct_summary.auth_asym_id",
                "_dssp_struct_summary.label_asym_id",
                "_dssp_struct_summary.asym_id",
            ],
            want_suffix=[".auth_asym_id", ".label_asym_id", ".asym_id"],
        )

        ss_idx = _find_tag_index(
            tags,
            want_exact=[
                "_dssp_struct_summary.secondary_structure",
                "_dssp_struct_summary.ss",
                "_dssp_struct_summary.structure",
            ],
            want_suffix=[".secondary_structure", ".ss", ".structure"],
        )

        beta: Dict[str, int] = {}
        if chain_idx is None or ss_idx is None:
            # category found but columns unexpected; keep going to fallback
            break

        for row in loop:
            chain = str(row[chain_idx]).strip()
            ss = str(row[ss_idx]).strip()
            if not chain or not ss:
                continue
            code = ss[0]
            if code in BETA_SS:
                beta[chain] = beta.get(chain, 0) + 1
        return beta

    # ---------- Fallback: broader heuristic (if primary fails) ----------
    beta: Dict[str, int] = {}
    for loop in block.loops:
        tags = [t.lower() for t in loop.tags]

        # try to find chain column
        chain_idx = None
        for i, t in enumerate(tags):
            if any(k in t for k in ("auth_asym_id", "label_asym_id", "asym_id", "chain_id")):
                chain_idx = i
                break

        # try to find ss column
        ss_idx = None
        for i, t in enumerate(tags):
            if any(k in t for k in ("secondary_structure", "sec_struct", "dssp", ".ss", " ss", "structure")):
                ss_idx = i
                break

        if chain_idx is None or ss_idx is None:
            continue

        try:
            for row in loop:
                chain = str(row[chain_idx]).strip()
                ss = str(row[ss_idx]).strip()
                if not chain or not ss:
                    continue
                code = ss[0]
                if code in BETA_SS:
                    beta[chain] = beta.get(chain, 0) + 1
            if beta:
                return beta
        except Exception:
            continue

    return beta


def process_one(cif_path: Path) -> tuple[str, bool, str, Dict[str, int], str]:
    """
    Returns: (pdb_id, success, mode, counts, err_msg)
    """
    pdb_id = cif_path.stem.upper()
    out_text = DSSP_DIR / f"{pdb_id}.dssp"
    out_cif = DSSP_DIR / f"{pdb_id}.dssp.cif"

    ok, out_used, mode, err = run_mkdssp_adaptive(cif_path, out_text, out_cif)
    if not ok:
        return pdb_id, False, mode, {}, err

    if mode == "dssp":
        counts = parse_dssp_text_beta_counts(out_used)
    else:
        counts = parse_mkdssp_mmcif_beta_counts(out_used)

    return pdb_id, True, mode, counts, ""


def main():
    cif_files = sorted(PDB_MODELS_DIR.glob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found in {PDB_MODELS_DIR}")
    if not MKDSSP_EXE.exists():
        raise FileNotFoundError(f"mkdssp.exe not found: {MKDSSP_EXE}")
    if not LIBCIFPP_DATA_DIR.exists():
        raise FileNotFoundError(f"LIBCIFPP_DATA_DIR not found: {LIBCIFPP_DATA_DIR}")

    print("===== STEP 5.1 (DSSP + beta counts) =====")
    print("CIF dir:", PDB_MODELS_DIR)
    print("DSSP dir:", DSSP_DIR)
    print("mkdssp:", MKDSSP_EXE)
    print("LIBCIFPP_DATA_DIR:", LIBCIFPP_DATA_DIR)
    print("Workers:", WORKERS)
    print("Timeout/file (sec):", MKDSSP_TIMEOUT_SEC)
    print("Total CIFs:", len(cif_files))

    # overwrite logs
    FAIL_LOG.write_text("", encoding="utf-8")
    ZEROCOUNT_LOG.write_text("", encoding="utf-8")

    # Write TSV header once, then append rows as files finish
    BETA_COUNTS_TSV.write_text("pdb_id\tchain_id\tbeta_residue_count\n", encoding="utf-8")
    tsv_lock = threading.Lock()

    ok = 0
    fail = 0
    zero = 0

    start_time = time.time()

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process_one, cif): cif for cif in cif_files}

        done = 0
        for fut in as_completed(futs):
            done += 1
            pdb_id, success, mode, counts, err = fut.result()

            if not success:
                fail += 1
                msg = f"[FAIL] {pdb_id} mode={mode}\n{err}\n{'-'*80}\n"
                print(f"[FAIL] {pdb_id}")
                with FAIL_LOG.open("a", encoding="utf-8") as f:
                    f.write(msg)
            else:
                ok += 1
                if not counts:
                    zero += 1
                    msg = f"[ZERO] {pdb_id} mode={mode} (no beta residues counted)\n"
                    with ZEROCOUNT_LOG.open("a", encoding="utf-8") as f:
                        f.write(msg)

                # append rows safely
                if counts:
                    with tsv_lock:
                        with BETA_COUNTS_TSV.open("a", encoding="utf-8", newline="\n") as out:
                            for chain, cnt in sorted(counts.items()):
                                out.write(f"{pdb_id}\t{chain}\t{cnt}\n")

            # frequent progress so it never looks “stuck”
            if done % 10 == 0 or done == len(cif_files):
                elapsed = time.time() - start_time
                print(
                    f"Progress: {done}/{len(cif_files)} "
                    f"(ok={ok}, fail={fail}, zero={zero}) "
                    f"elapsed={elapsed:.1f}s"
                )

    print("\n===== STEP 5.1 DONE =====")
    print("Total CIFs:", len(cif_files))
    print("DSSP OK:", ok)
    print("DSSP FAIL:", fail)
    print("DSSP ZERO beta-count files:", zero)
    print("Wrote:", BETA_COUNTS_TSV)
    print("Fail log:", FAIL_LOG)
    print("Zero-count log:", ZEROCOUNT_LOG)


if __name__ == "__main__":
    main()

# Run:
#   .\venv\Scripts\python.exe -m pip install gemmi
#   .\venv\Scripts\python.exe .\step5_1_make_beta_counts_dssp.py
