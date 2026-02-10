from __future__ import annotations
from pathlib import Path
import shutil
import subprocess

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "outputs"

# Input from your project (may contain spaces — that's fine)
FASTA_IN = OUT_DIR / "all_chains_beta_only.fasta"

# Scratch workdir with NO spaces
SCRATCH = Path(r"C:\mmseqs_beta_only")
SCRATCH.mkdir(parents=True, exist_ok=True)

FASTA = SCRATCH / "all_chains_beta_only.fasta"
WORK = SCRATCH / "work"
WORK.mkdir(exist_ok=True)

# New output (won't overwrite old files)
OUT_TSV = OUT_DIR / "mmseqs_clusters_beta_only.tsv"

MMSEQS = "mmseqs"   # if not in PATH, use full path to mmseqs.exe
THREADS = 8

DB = WORK / "chainsDB"
TMP = WORK / "tmp"
CLUST = WORK / "clustRes"

def run(cmd: list[str]):
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    if not FASTA_IN.exists():
        raise FileNotFoundError(f"Missing FASTA: {FASTA_IN}")

    # Safety: don't overwrite TSV
    if OUT_TSV.exists():
        raise FileExistsError(f"{OUT_TSV} already exists. Rename/delete it to rerun safely.")

    # Copy FASTA into no-space path
    shutil.copy2(FASTA_IN, FASTA)

    # createdb (force threads to 8)
    run([MMSEQS, "createdb", str(FASTA), str(DB), "--threads", str(THREADS)])

    # cluster
    run([
        MMSEQS, "cluster",
        str(DB), str(CLUST), str(TMP),
        "--min-seq-id", "0.9",
        "-c", "0.8",
        "--cov-mode", "1",
        "--threads", str(THREADS)
    ])

    # export rep-member TSV
    run([MMSEQS, "createtsv", str(DB), str(DB), str(CLUST), str(OUT_TSV)])

    print("\nDONE")
    print("Wrote:", OUT_TSV)
    print("Scratch workdir:", WORK)

if __name__ == "__main__":
    main()
# .\venv\Scripts\python.exe .\step_mmseqs_cluster_beta_only_nospace.py
