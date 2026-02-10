import time, random
from pathlib import Path
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices

FASTA = Path(r"C:\Users\fd02629\Desktop\Protein Experiment\outputs\all_chains.fasta")
GAP_OPEN, GAP_EXTEND = -10.0, -0.5
blosum62 = substitution_matrices.load("BLOSUM62")

records = list(SeqIO.parse(str(FASTA), "fasta"))
seqs = [str(r.seq) for r in records]

# pick random 100 pairs
pairs = []
n = len(seqs)
for _ in range(100):
    i = random.randrange(n)
    j = random.randrange(n)
    while j == i:
        j = random.randrange(n)
    pairs.append((seqs[i], seqs[j]))

t0 = time.time()
for a, b in pairs:
    pairwise2.align.globalds(a, b, blosum62, GAP_OPEN, GAP_EXTEND, one_alignment_only=True)
t1 = time.time()

total = t1 - t0
print("100 alignments seconds:", total)
print("seconds per alignment:", total/100)
