from Bio.Align import substitution_matrices
from Bio import pairwise2

blosum62 = substitution_matrices.load("BLOSUM62")

def nw_identity_coverage(seq1, seq2, gap_open=-10.0, gap_extend=-0.5):
    aln = pairwise2.align.globalds(seq1, seq2, blosum62, gap_open, gap_extend, one_alignment_only=True)[0]
    a1, a2, score, start, end = aln
    matches = 0
    aligned_res = 0
    for c1, c2 in zip(a1, a2):
        if c1 != "-" and c2 != "-":
            aligned_res += 1
            if c1 == c2:
                matches += 1
    identity = matches / aligned_res if aligned_res else 0.0
    coverage = aligned_res / min(len(seq1), len(seq2)) if min(len(seq1), len(seq2)) else 0.0
    return identity, coverage

s1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQLR"
s2 = "MKTAYIAKQRQISFVKSHFSRQLAERLGLIEVQLR"

iden, cov = nw_identity_coverage(s1, s2)
print("identity:", iden)
print("coverage:", cov)
print("similar?", iden >= 0.90 and cov >= 0.80)
