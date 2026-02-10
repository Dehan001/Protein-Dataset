from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple
import subprocess
import re

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"          # your CIF files
OUT_DIR = BASE / "outputs"
DSSP_DIR = BASE / "dssp_out"

# IMPORTANT: point to mkdssp.exe you used (working)
MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")

# Output TSV
BETA_COUNTS_TSV = OUT_DIR / "beta_counts.tsv"
# ====================================================

OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

# DSSP lines after header look like:
#  #  RESIDUE AA STRUCTURE ...
# We parse: chain ID column and structure column.
# DSSP "E" = beta strand, "B" = beta bridge.
BETA_SS = {"E", "B"}


def run_dssp(cif_path: Path, dssp_path: Path) -> bool:
    """Run mkdssp on one cif -> dssp."""
    try:
        # mkdssp input output are positional args in your build
        subprocess.run(
            [str(MKDSSP_EXE), str(cif_path), str(dssp_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"[FAIL] mkdssp: {cif_path.name}")
        # Print only short error tail
        err = (e.stderr or "").strip().splitlines()[-5:]
        for line in err:
            print("   ", line)
        return False


def parse_dssp_beta_counts(dssp_path: Path) -> Dict[str, int]:
    """
    Return dict: chain_id -> beta_residue_count
    (counts residues assigned to E or B)
    """
    beta: Dict[str, int] = {}
    with dssp_path.open("r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # DSSP header ends at line starting with "  #  RESIDUE"
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

        # DSSP fixed-width:
        # chain id at position 11 (0-based index 11) in classic DSSP output.
        # secondary structure at position 16.
        # These positions are consistent for mkdssp outputs.
        chain = line[11].strip()
        ss = line[16].strip()

        if not chain:
            continue
        if ss in BETA_SS:
            beta[chain] = beta.get(chain, 0) + 1

    return beta


def main():
    cif_files = sorted(PDB_MODELS_DIR.glob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found in {PDB_MODELS_DIR}")

    if not MKDSSP_EXE.exists():
        raise FileNotFoundError(f"mkdssp.exe not found: {MKDSSP_EXE}")

    total = 0
    ok = 0
    fail = 0

    # aggregate: (pdb_id, chain) -> beta_count
    all_counts: Dict[Tuple[str, str], int] = {}

    for cif in cif_files:
        total += 1
        pdb_id = cif.stem.upper()
        dssp_path = DSSP_DIR / f"{pdb_id}.dssp"

        if not dssp_path.exists():
            success = run_dssp(cif, dssp_path)
            if not success:
                fail += 1
                continue

        beta_by_chain = parse_dssp_beta_counts(dssp_path)
        ok += 1

        for chain, cnt in beta_by_chain.items():
            all_counts[(pdb_id, chain)] = cnt

        if total % 50 == 0:
            print(f"Processed {total}/{len(cif_files)} (ok={ok}, fail={fail})")

    # write beta_counts.tsv
    with BETA_COUNTS_TSV.open("w", encoding="utf-8", newline="\n") as out:
        out.write("pdb_id\tchain_id\tbeta_residue_count\n")
        for (pdb_id, chain), cnt in sorted(all_counts.items()):
            out.write(f"{pdb_id}\t{chain}\t{cnt}\n")

    print("\n===== STEP 5.1 DONE =====")
    print("Total CIFs:", len(cif_files))
    print("DSSP OK:", ok)
    print("DSSP FAIL:", fail)
    print("Wrote:", BETA_COUNTS_TSV)


if __name__ == "__main__":
    main()
