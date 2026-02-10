from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import os

import gemmi  # pip install gemmi

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"         # your *.cif coordinate files
DSSP_DIR = BASE / "dssp_out"                 # where we store mkdssp outputs (for records)
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

# Working mkdssp.exe
MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")

# Where your dictionary lives (you found it already)
LIBCIFPP_DATA_DIR = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\share\libcifpp")

# Threads
WORKERS = 8

# Output
BETA_COUNTS_TSV = OUT_DIR / "beta_counts.tsv"

# DSSP beta codes
BETA_SS = {"E", "B"}
# ====================================================


def run_cmd(cmd: list[str], env: dict) -> tuple[int, str, str]:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    return p.returncode, p.stdout, p.stderr


def run_mkdssp_adaptive(cif_path: Path, out_text_path: Path, out_cif_path: Path) -> tuple[bool, Path, str]:
    """
    Try to run mkdssp producing classic DSSP text first.
    If mkdssp complains about old DSSP format, rerun to produce mmCIF output.
    Returns: (ok, output_path_used, mode) where mode is "dssp" or "mmcif".
    """

    env = os.environ.copy()
    env["LIBCIFPP_DATA_DIR"] = str(LIBCIFPP_DATA_DIR)

    # 1) Try classic DSSP text output (positional args: input output)
    rc, so, se = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_text_path)], env)

    if rc == 0 and out_text_path.exists():
        return True, out_text_path, "dssp"

    # If it fails due to old DSSP format, switch to mmCIF output
    if "old DSSP format" in (se or "") or "mmCIF format instead" in (se or ""):
        # Many mkdssp builds support --output-format; try it.
        rc2, so2, se2 = run_cmd(
            [str(MKDSSP_EXE), "--output-format=mmcif", str(cif_path), str(out_cif_path)],
            env,
        )
        if rc2 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif"

        # Fallback: some builds auto-detect from extension
        rc3, so3, se3 = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_cif_path)], env)
        if rc3 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif"

        return False, out_cif_path, "mmcif"

    # Other failures
    return False, out_text_path, "dssp"


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


def parse_mkdssp_mmcif_beta_counts(mmcif_path: Path) -> Dict[str, int]:
    """
    Heuristic mmCIF parser:
    Find a loop that contains:
      - a chain identifier column (auth_asym_id / label_asym_id / asym_id)
      - a secondary structure code column containing one-letter DSSP codes
    Then count E/B per chain.
    """
    doc = gemmi.cif.read_file(str(mmcif_path))
    block = doc.sole_block()

    beta: Dict[str, int] = {}

    for loop in block.loops:
        tags = [t.lower() for t in loop.tags]

        # chain column candidates
        chain_idx = None
        for key in ("auth_asym_id", "label_asym_id", "asym_id", "chain_id"):
            for i, t in enumerate(tags):
                if t.endswith(key) or key in t:
                    chain_idx = i
                    break
            if chain_idx is not None:
                break

        # ss column candidates
        ss_idx = None
        # common-ish names across DSSP mmCIF exports
        for key in ("sec_struct", "secondary_structure", "ss", "structure", "dssp"):
            for i, t in enumerate(tags):
                if any(k in t for k in ("sec", "secondary", "ss", "struct", "dssp")):
                    # prefer something that looks like it stores the one-letter code
                    if key in t:
                        ss_idx = i
                        break
            if ss_idx is not None:
                break

        if chain_idx is None or ss_idx is None:
            continue

        # try reading rows
        try:
            for row in loop:
                chain = str(row[chain_idx]).strip()
                ss = str(row[ss_idx]).strip()
                if not chain or not ss:
                    continue
                # usually one-letter; if longer, take first char
                code = ss[0]
                if code in BETA_SS:
                    beta[chain] = beta.get(chain, 0) + 1
            # if we successfully counted something, stop early
            if beta:
                return beta
        except Exception:
            continue

    return beta


def process_one(cif_path: Path) -> tuple[str, bool, int, Dict[str, int]]:
    pdb_id = cif_path.stem.upper()

    out_text = DSSP_DIR / f"{pdb_id}.dssp"
    out_cif = DSSP_DIR / f"{pdb_id}.dssp.cif"

    ok, out_used, mode = run_mkdssp_adaptive(cif_path, out_text, out_cif)
    if not ok:
        return pdb_id, False, 0, {}

    if mode == "dssp":
        counts = parse_dssp_text_beta_counts(out_used)
    else:
        counts = parse_mkdssp_mmcif_beta_counts(out_used)

    return pdb_id, True, len(counts), counts


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
    print("Total CIFs:", len(cif_files))

    all_counts: Dict[Tuple[str, str], int] = {}
    ok = 0
    fail = 0

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process_one, cif): cif for cif in cif_files}

        done = 0
        for fut in as_completed(futs):
            done += 1
            pdb_id, success, n_chains, counts = fut.result()

            if not success:
                fail += 1
                print(f"[FAIL] {pdb_id}")
            else:
                ok += 1
                for chain, cnt in counts.items():
                    all_counts[(pdb_id, chain)] = cnt

            if done % 50 == 0:
                print(f"Progress: {done}/{len(cif_files)} (ok={ok}, fail={fail})")

    with BETA_COUNTS_TSV.open("w", encoding="utf-8", newline="\n") as out:
        out.write("pdb_id\tchain_id\tbeta_residue_count\n")
        for (pdb_id, chain), cnt in sorted(all_counts.items()):
            out.write(f"{pdb_id}\t{chain}\t{cnt}\n")

    print("\n===== STEP 5.1 DONE =====")
    print("Total CIFs:", len(cif_files))
    print("DSSP OK:", ok)
    print("DSSP FAIL:", fail)
    print("Beta rows written:", len(all_counts))
    print("Wrote:", BETA_COUNTS_TSV)


if __name__ == "__main__":
    main()


# .\venv\Scripts\python.exe -m pip install gemmi
# .\venv\Scripts\python.exe .\step5_1_beta_counts_dssp.py
