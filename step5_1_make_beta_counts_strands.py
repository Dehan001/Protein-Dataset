from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import os
import threading
import time

import gemmi  # pip install gemmi
from Bio import SeqIO  # pip install biopython

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"         # your *.cif coordinate files
DSSP_DIR = BASE / "dssp_out"                 # where we store mkdssp outputs
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

# FASTA with IDs like: PDBID_EMD-xxxxx_A
ALL_CHAINS_FASTA = OUT_DIR / "all_chains.fasta"

# mkdssp that you confirmed works when called by full path
MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")

# Important for the "mmcif_pdbx.dic" issue
LIBCIFPP_DATA_DIR = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\share\libcifpp")

WORKERS = 8
MKDSSP_TIMEOUT_SEC = 900  # 15 min per file

BETA_COUNTS_TSV = OUT_DIR / "beta_counts.tsv"
FAIL_LOG = OUT_DIR / "dssp_failures.log"
ZEROCOUNT_LOG = OUT_DIR / "dssp_zero_counts.log"

# DSSP beta codes
BETA_SS = {"E", "B"}
# ====================================================


def build_fasta_id_map(fasta_path: Path) -> Dict[Tuple[str, str], list[str]]:
    """
    Map (PDB_ID, CHAIN) -> [FASTA_IDs]
    FASTA IDs expected: PDBID_EMD-xxxxx_CHAIN
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing {fasta_path}")

    m: Dict[Tuple[str, str], list[str]] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        fid = rec.id.strip()
        parts = fid.split("_")
        if len(parts) < 2:
            continue
        pdb = parts[0].upper()
        chain = parts[-1]  # last token is chain ID
        key = (pdb, chain)
        m.setdefault(key, []).append(fid)
    return m


def run_cmd(cmd: list[str], env: dict, timeout_sec: int) -> tuple[int, str, str, bool]:
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
    Try classic DSSP text first; if mkdssp says "use mmCIF format", rerun to produce mmCIF output.
    Returns: (ok, output_path_used, mode, err_msg) where mode in {"dssp","mmcif"}
    """
    env = os.environ.copy()
    env["LIBCIFPP_DATA_DIR"] = str(LIBCIFPP_DATA_DIR)

    # 1) classic DSSP text: mkdssp input output
    rc, so, se, _ = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_text_path)], env, MKDSSP_TIMEOUT_SEC)
    if rc == 0 and out_text_path.exists():
        return True, out_text_path, "dssp", ""

    err = (se or "") + (so or "")

    # 2) switch to mmCIF output if needed
    if "old DSSP format" in err or "mmCIF format instead" in err:
        rc2, so2, se2, _ = run_cmd(
            [str(MKDSSP_EXE), "--output-format=mmcif", str(cif_path), str(out_cif_path)],
            env,
            MKDSSP_TIMEOUT_SEC,
        )
        if rc2 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif", ""

        # fallback: infer from extension
        rc3, so3, se3, _ = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_cif_path)], env, MKDSSP_TIMEOUT_SEC)
        if rc3 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif", ""

        return False, out_cif_path, "mmcif", (se2 or "") + (so2 or "") + (se3 or "") + (so3 or "")

    return False, out_text_path, "dssp", err


def count_beta_strands_from_dssp_text(dssp_path: Path) -> Dict[str, int]:
    """
    Count β-STRANDS (segments), not residues.
    For each chain: increment when SS enters {E,B} from non-{E,B}.
    """
    lines = dssp_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    start = None
    for idx, line in enumerate(lines):
        if line.startswith("  #  RESIDUE"):
            start = idx + 1
            break
    if start is None:
        return {}

    beta_strands: Dict[str, int] = {}
    prev_is_beta: Dict[str, bool] = {}

    for line in lines[start:]:
        if len(line) < 17:
            continue
        chain = line[11].strip()
        ss = line[16].strip()
        if not chain:
            continue

        is_beta = (ss in BETA_SS)
        was_beta = prev_is_beta.get(chain, False)
        if is_beta and not was_beta:
            beta_strands[chain] = beta_strands.get(chain, 0) + 1
        prev_is_beta[chain] = is_beta

    return beta_strands


def parse_mkdssp_mmcif_beta_strands(mmcif_path: Path) -> Dict[str, int]:
    """
    Robust mmCIF parser for mkdssp output:
    Count β-strand segments per chain (transitions into E/B).
    """
    doc = gemmi.cif.read_file(str(mmcif_path))
    block = doc.sole_block()

    best: Dict[str, int] = {}

    # Try loops; look for DSSP-like categories first
    for loop in block.loops:
        tags = [t.lower() for t in loop.tags]

        # heuristic: if it smells like DSSP output
        if not any(("dssp" in t or "secondary_structure" in t or "sec_struct" in t) for t in tags):
            continue

        # find chain-ish column
        chain_idx = None
        for i, t in enumerate(tags):
            if any(k in t for k in ("auth_asym_id", "label_asym_id", "asym_id", "chain_id")):
                chain_idx = i
                break

        # find ss-ish column
        ss_idx = None
        for i, t in enumerate(tags):
            if any(k in t for k in ("secondary_structure", "sec_struct", "structure", ".ss", "dssp")):
                ss_idx = i
                break

        if chain_idx is None or ss_idx is None:
            continue

        # count transitions
        beta_strands: Dict[str, int] = {}
        prev_is_beta: Dict[str, bool] = {}

        try:
            for row in loop:
                chain = str(row[chain_idx]).strip()
                ss = str(row[ss_idx]).strip()
                if not chain or not ss:
                    continue
                code = ss[0]
                is_beta = (code in BETA_SS)
                was_beta = prev_is_beta.get(chain, False)
                if is_beta and not was_beta:
                    beta_strands[chain] = beta_strands.get(chain, 0) + 1
                prev_is_beta[chain] = is_beta
        except Exception:
            continue

        if beta_strands:
            return beta_strands

        best = beta_strands

    return best


def process_one(cif_path: Path) -> tuple[str, bool, str, Dict[str, int], str]:
    pdb_id = cif_path.stem.upper()
    out_text = DSSP_DIR / f"{pdb_id}.dssp"
    out_cif = DSSP_DIR / f"{pdb_id}.dssp.cif"

    ok, out_used, mode, err = run_mkdssp_adaptive(cif_path, out_text, out_cif)
    if not ok:
        return pdb_id, False, mode, {}, err

    if mode == "dssp":
        counts = count_beta_strands_from_dssp_text(out_used)
    else:
        counts = parse_mkdssp_mmcif_beta_strands(out_used)

    return pdb_id, True, mode, counts, ""


def main():
    cif_files = sorted(PDB_MODELS_DIR.glob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found in {PDB_MODELS_DIR}")
    if not MKDSSP_EXE.exists():
        raise FileNotFoundError(f"mkdssp.exe not found: {MKDSSP_EXE}")
    if not LIBCIFPP_DATA_DIR.exists():
        raise FileNotFoundError(f"LIBCIFPP_DATA_DIR not found: {LIBCIFPP_DATA_DIR}")

    fasta_map = build_fasta_id_map(ALL_CHAINS_FASTA)

    print("===== STEP 5.1 (DSSP -> beta STRANDS per FASTA chain) =====")
    print("Total CIFs:", len(cif_files))
    print("Workers:", WORKERS)

    FAIL_LOG.write_text("", encoding="utf-8")
    ZEROCOUNT_LOG.write_text("", encoding="utf-8")
    BETA_COUNTS_TSV.write_text("fasta_id\tbeta_strands\n", encoding="utf-8")
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
            pdb_id, success, mode, counts_by_chain, err = fut.result()

            if not success:
                fail += 1
                with FAIL_LOG.open("a", encoding="utf-8") as f:
                    f.write(f"[FAIL] {pdb_id} mode={mode}\n{err}\n{'-'*80}\n")
            else:
                ok += 1
                if not counts_by_chain:
                    zero += 1
                    with ZEROCOUNT_LOG.open("a", encoding="utf-8") as f:
                        f.write(f"[ZERO] {pdb_id} mode={mode} (no beta strands found)\n")
                else:
                    # write per FASTA chain ID (PDB+CHAIN -> possibly multiple FASTA IDs)
                    rows = []
                    for chain, b in counts_by_chain.items():
                        key = (pdb_id, chain)
                        for fasta_id in fasta_map.get(key, []):
                            rows.append((fasta_id, b))

                    if rows:
                        with tsv_lock:
                            with BETA_COUNTS_TSV.open("a", encoding="utf-8", newline="\n") as out:
                                for fasta_id, b in sorted(rows):
                                    out.write(f"{fasta_id}\t{b}\n")

            if done % 10 == 0 or done == len(cif_files):
                elapsed = time.time() - start_time
                print(f"Progress: {done}/{len(cif_files)} (ok={ok}, fail={fail}, zero={zero}) elapsed={elapsed:.1f}s")

    print("\n===== STEP 5.1 DONE =====")
    print("OK:", ok, "FAIL:", fail, "ZERO:", zero)
    print("Wrote:", BETA_COUNTS_TSV)
    print("Fail log:", FAIL_LOG)
    print("Zero log:", ZEROCOUNT_LOG)


if __name__ == "__main__":
    main()


# in your project folder
# .\venv\Scripts\python.exe -m pip install gemmi biopython
# .\venv\Scripts\python.exe .\step5_1_make_beta_counts_strands.py
