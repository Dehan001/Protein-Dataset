from __future__ import annotations
from pathlib import Path
import subprocess
import os

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"

FASTA = OUT / "all_chains_ss_nonzero.fasta"

# IMPORTANT: a work directory with NO SPACES
WORK = Path(r"C:\mmseqs_tmp\ss_nonzero")     # <-- change if you want, but keep no spaces
DB = WORK / "chainsDB"
CLUST = WORK / "clustRes"
TMP = WORK / "tmp"

OUT_TSV = OUT / "mmseqs_clusters_ss_nonzero.tsv"

MMSEQS = "mmseqs"
THREADS = 8


def run(cmd: list[str]):
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    if not FASTA.exists():
        raise FileNotFoundError(FASTA)

    WORK.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)

    # (optional) clean old run outputs but DO NOT touch your outputs folder
    # If you want totally fresh run each time:
    # for p in [DB, CLUST]:
    #     if p.exists():
    #         # mmseqs creates multiple files with suffixes; remove by prefix
    #         for f in WORK.glob(p.name + "*"):
    #             try: f.unlink()
    #             except: pass

    run([MMSEQS, "createdb", str(FASTA), str(DB)])

    run([
        MMSEQS, "cluster",
        str(DB), str(CLUST), str(TMP),
        "--min-seq-id", "0.9",
        "-c", "0.8",
        "--cov-mode", "1",
        "--threads", str(THREADS),
    ])

    run([MMSEQS, "createtsv", str(DB), str(DB), str(CLUST), str(OUT_TSV)])

    print("Wrote:", OUT_TSV)
    print("Workdir:", WORK)


if __name__ == "__main__":
    main()

# .\venv\Scripts\python.exe .\step2_mmseqs_cluster_ss_nonzero.py
