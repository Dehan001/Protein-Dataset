from pathlib import Path
from Bio import SeqIO

# ===================== PATHS =====================
BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "outputs"

IN_FASTA = OUT_DIR / "all_chains.fasta"
BETA_TSV = OUT_DIR / "beta_counts.tsv"

OUT_FASTA = OUT_DIR / "all_chains_beta_only.fasta"
REMOVED_TXT = OUT_DIR / "chains_removed_no_beta.txt"

# ===================== LOAD BETA COUNTS =====================
beta_counts = {}

with BETA_TSV.open("r", encoding="utf-8") as f:
    header = f.readline().strip().split("\t")
    if header != ["fasta_id", "beta_strands"]:
        raise ValueError("beta_counts.tsv must have columns: fasta_id<TAB>beta_strands")

    for line in f:
        fid, b = line.rstrip("\n").split("\t")
        beta_counts[fid] = int(b)

# ===================== FILTER FASTA =====================
kept = []
removed = []

for rec in SeqIO.parse(str(IN_FASTA), "fasta"):
    fid = rec.id
    b = beta_counts.get(fid, 0)

    if b > 0:
        kept.append(rec)
    else:
        removed.append((fid, b))

# ===================== WRITE OUTPUTS =====================
SeqIO.write(kept, str(OUT_FASTA), "fasta")

with REMOVED_TXT.open("w", encoding="utf-8", newline="\n") as f:
    f.write("fasta_id\tbeta_strands\n")
    for fid, b in removed:
        f.write(f"{fid}\t{b}\n")

# ===================== SUMMARY =====================
print("===== REMOVE NO-BETA CHAINS =====")
print(f"Input FASTA chains      : {len(kept) + len(removed)}")
print(f"Kept (beta > 0)         : {len(kept)}")
print(f"Removed (beta == 0)     : {len(removed)}")
print("Wrote:")
print(" -", OUT_FASTA)
print(" -", REMOVED_TXT)

# .\venv\Scripts\python.exe .\step_remove_no_beta_chains.py
