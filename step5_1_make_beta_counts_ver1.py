from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, List, Optional
import subprocess
import multiprocessing as mp

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"          # your CIF files
OUT_DIR = BASE / "outputs"
DSSP_DIR = BASE / "dssp_out"

# IMPORTANT: point to mkdssp.exe you used (working)
MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")

# Output TSV
BETA_COUNTS_TSV = OUT_DIR / "beta_counts.tsv"

# Parallelism
N_WORKERS = 8          # as requested
CHUNKSIZE = 10         # how many CIFs each worker grabs per scheduling chunk

# DSSP "E" = beta strand, "B" = beta bridge.
BETA_SS = {"E", "B"}
# ====================================================

OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

# ---------- Worker Globals ----------
G_MKDSSP: Optional[str] = None
G_DSSP_DIR: Optional[str] = None
# -----------------------------------


def _init_worker(mkdssp_exe: str, dssp_dir: str):
    global G_MKDSSP, G_DSSP_DIR
    G_MKDSSP = mkdssp_exe
    G_DSSP_DIR = dssp_dir


def run_dssp(cif_path: Path, dssp_path: Path) -> Tuple[bool, str]:
    """Run mkdssp on one cif -> dssp. Returns (ok, error_tail)."""
    try:
        subprocess.run(
            [str(G_MKDSSP), str(cif_path), str(dssp_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return True, ""
    except subprocess.CalledProcessError as e:
        err_lines = (e.stderr or "").strip().splitlines()
        tail = "\n".join(err_lines[-5:]) if err_lines else "(no stderr)"
        return False, tail
    except Exception as e:
        return False, str(e)


def parse_dssp_beta_counts(dssp_path: Path) -> Dict[str, int]:
    """
    Return dict: chain_id -> beta_residue_count
    (counts residues assigned to E or B)
    """
    beta: Dict[str, int] = {}

    with dssp_path.open("r", encoding="utf-8", errors="ignore") as f:
        in_data = False
        for line in f:
            if line.startswith("  #  RESIDUE"):
                in_data = True
                continue
            if not in_data:
                continue
            if len(line) < 17:
                continue

            chain = line[11].strip()
            ss = line[16].strip()

            if not chain:
                continue
            if ss in BETA_SS:
                beta[chain] = beta.get(chain, 0) + 1

    return beta


def process_one_cif(cif_path_str: str) -> Tuple[str, bool, str, Dict[Tuple[str, str], int]]:
    """
    Worker task:
    - ensure DSSP exists (create if missing)
    - parse DSSP and return {(pdb_id, chain): beta_count}
    Returns:
      pdb_id, ok, error_message, counts_dict
    """
    cif = Path(cif_path_str)
    pdb_id = cif.stem.upper()

    dssp_path = Path(G_DSSP_DIR) / f"{pdb_id}.dssp"

    # create DSSP if missing
    if not dssp_path.exists():
        ok, err = run_dssp(cif, dssp_path)
        if not ok:
            return pdb_id, False, err, {}

    # parse DSSP
    try:
        beta_by_chain = parse_dssp_beta_counts(dssp_path)
    except Exception as e:
        return pdb_id, False, f"parse error: {e}", {}

    out: Dict[Tuple[str, str], int] = {}
    for chain, cnt in beta_by_chain.items():
        out[(pdb_id, chain)] = cnt

    return pdb_id, True, "", out


def main():
    cif_files = sorted(PDB_MODELS_DIR.glob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found in {PDB_MODELS_DIR}")

    if not MKDSSP_EXE.exists():
        raise FileNotFoundError(f"mkdssp.exe not found: {MKDSSP_EXE}")

    print("===== STEP 5.1 (DSSP + beta counts) =====")
    print("CIF dir:", PDB_MODELS_DIR)
    print("DSSP dir:", DSSP_DIR)
    print("mkdssp:", MKDSSP_EXE)
    print("Workers:", N_WORKERS)
    print("Total CIFs:", len(cif_files))

    all_counts: Dict[Tuple[str, str], int] = {}
    ok = 0
    fail = 0

    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=N_WORKERS,
        initializer=_init_worker,
        initargs=(str(MKDSSP_EXE), str(DSSP_DIR)),
        maxtasksperchild=200,   # helps stability for long runs on Windows
    ) as pool:

        for idx, (pdb_id, success, err, counts) in enumerate(
            pool.imap_unordered(process_one_cif, [str(p) for p in cif_files], chunksize=CHUNKSIZE),
            start=1
        ):
            if success:
                ok += 1
                all_counts.update(counts)
            else:
                fail += 1
                print(f"[FAIL] {pdb_id}")
                if err:
                    for line in err.splitlines()[-5:]:
                        print("   ", line)

            if idx % 50 == 0:
                print(f"Progress: {idx}/{len(cif_files)} (ok={ok}, fail={fail})")

    # write beta_counts.tsv
    with BETA_COUNTS_TSV.open("w", encoding="utf-8", newline="\n") as out:
        out.write("pdb_id\tchain_id\tbeta_residue_count\n")
        for (pdb_id, chain), cnt in sorted(all_counts.items()):
            out.write(f"{pdb_id}\t{chain}\t{cnt}\n")

    print("\n===== DONE =====")
    print("Total CIFs:", len(cif_files))
    print("DSSP OK:", ok)
    print("DSSP FAIL:", fail)
    print("Beta rows written:", len(all_counts))
    print("Wrote:", BETA_COUNTS_TSV)


if __name__ == "__main__":
    main()
