from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple
from Bio import SeqIO

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"

# INPUTS
IN_FASTA = OUT / "all_chains.fasta"  # or your clustering-before-annotation FASTA
SS_COUNTS = OUT / "ss_counts_20260203_125637.tsv"  # <-- change to your real file

# OUTPUTS (new, no overwrite)
TAG = "ss_nonzero"
OUT_FASTA = OUT / f"all_chains_{TAG}.fasta"
OUT_REMOVED = OUT / f"chains_removed_all_ss_zero_{TAG}.tsv"
OUT_SUMMARY = OUT / f"filter_summary_{TAG}.txt"


def load_ss(path: Path) -> Dict[str, Tuple[int, int, int, int]]:
    ss: Dict[str, Tuple[int, int, int, int]] = {}
    with path.open("r", encoding="utf-8") as f:
        _ = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            fid = parts[0].strip()
            try:
                br = int(parts[1])
                hr = int(parts[2])
                bs = int(parts[3])
                hx = int(parts[4])
            except ValueError:
                br = hr = bs = hx = 0
            ss[fid] = (br, hr, bs, hx)
    return ss


def main():
    if not IN_FASTA.exists():
        raise FileNotFoundError(IN_FASTA)
    if not SS_COUNTS.exists():
        raise FileNotFoundError(SS_COUNTS)

    ss = load_ss(SS_COUNTS)

    kept = []
    removed = []

    total = 0
    for rec in SeqIO.parse(str(IN_FASTA), "fasta"):
        total += 1
        fid = rec.id
        br, hr, bs, hx = ss.get(fid, (0, 0, 0, 0))

        # Remove only if EVERYTHING is zero
        if (br == 0 and hr == 0 and bs == 0 and hx == 0):
            removed.append((fid, br, hr, bs, hx))
        else:
            kept.append(rec)

    SeqIO.write(kept, str(OUT_FASTA), "fasta")

    with OUT_REMOVED.open("w", encoding="utf-8", newline="\n") as f:
        f.write("fasta_id\tbeta_res\thelix_res\tbeta_strands\thelices\n")
        for fid, br, hr, bs, hx in removed:
            f.write(f"{fid}\t{br}\t{hr}\t{bs}\t{hx}\n")

    OUT_SUMMARY.write_text(
        f"Input FASTA chains: {total}\n"
        f"Kept (any SS > 0): {len(kept)}\n"
        f"Removed (all SS == 0): {len(removed)}\n"
        f"Zero-rate: {len(removed)/max(total,1):.4f}\n"
        f"SS source: {SS_COUNTS.name}\n"
        f"Wrote: {OUT_FASTA.name}\n"
        f"Wrote: {OUT_REMOVED.name}\n",
        encoding="utf-8"
    )

    print("Wrote:", OUT_FASTA)
    print("Wrote:", OUT_REMOVED)
    print("Wrote:", OUT_SUMMARY)


if __name__ == "__main__":
    main()

# .\venv\Scripts\python.exe .\step1_filter_zero_ss.py
