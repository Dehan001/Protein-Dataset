from pathlib import Path
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices

# --------- SETTINGS (match your instruction) ----------
FASTA_PATH = Path(r"C:\Protein\outputs\all_chains.fasta")
OUT_DIR = Path(r"C:\Protein\outputs")

MIN_IDENTITY = 0.90
MIN_COVERAGE = 0.80

# Pairwise2 globalds needs gap penalties. These are common defaults.
GAP_OPEN = -10.0
GAP_EXTEND = -0.5

# Speed filter: skip if max_len/min_len > LENGTH_RATIO
# If you want "no filters at all", set this very high (e.g., 1000).
# LENGTH_RATIO = 1.25
LENGTH_RATIO = 1.25

# Optional: stop early for testing; set to None to run all
LIMIT_N = 50
# ------------------------------------------------------

OUT_DIR.mkdir(parents=True, exist_ok=True)

blosum62 = substitution_matrices.load("BLOSUM62")


class UnionFind:
    """
    Union-Find that enforces: higher ID attaches to lower ID (lower ID becomes rep).
    """
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union_lower_wins(self, a: int, b: int) -> bool:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return False
        # lower root wins
        if ra < rb:
            self.parent[rb] = ra
        else:
            self.parent[ra] = rb
        return True


def nw_identity_coverage(seq1: str, seq2: str):
    """
    Needleman-Wunsch global alignment using BLOSUM62 via pairwise2.globalds
    identity = matches / aligned_residues
    coverage = aligned_residues / min(len(seq1), len(seq2))
    """
    alns = pairwise2.align.globalds(
        seq1, seq2, blosum62, GAP_OPEN, GAP_EXTEND, one_alignment_only=True
    )
    if not alns:
        return 0.0, 0.0

    a1, a2, score, start, end = alns[0]

    matches = 0
    aligned = 0
    for c1, c2 in zip(a1, a2):
        if c1 != "-" and c2 != "-":
            aligned += 1
            if c1 == c2:
                matches += 1

    if aligned == 0:
        return 0.0, 0.0

    identity = matches / aligned
    coverage = aligned / min(len(seq1), len(seq2))
    return identity, coverage


def main():
    records = list(SeqIO.parse(str(FASTA_PATH), "fasta"))
    if LIMIT_N is not None:
        records = records[:LIMIT_N]

    ids = [r.id for r in records]
    seqs = [str(r.seq) for r in records]
    lens = [len(s) for s in seqs]
    n = len(seqs)

    print(f"Loaded sequences: {n}")
    print("Running full all-vs-all NW clustering with:")
    print(f"  MIN_IDENTITY={MIN_IDENTITY}, MIN_COVERAGE={MIN_COVERAGE}")
    print(f"  GAP_OPEN={GAP_OPEN}, GAP_EXTEND={GAP_EXTEND}")
    print(f"  LENGTH_RATIO filter={LENGTH_RATIO}")

    uf = UnionFind(n)

    # Progress counters
    tested = 0
    merged = 0

    # Main all-vs-all loop
    for i in range(n):
        s1 = seqs[i]
        l1 = lens[i]

        if i % 25 == 0:
            print(f"i={i}/{n} | tested={tested} | merged={merged}")

        for j in range(i + 1, n):
            l2 = lens[j]

            # quick length ratio filter (safe for 80% coverage threshold)
            mn = min(l1, l2)
            mx = max(l1, l2)
            if mx / mn > LENGTH_RATIO:
                continue

            tested += 1
            identity, coverage = nw_identity_coverage(s1, seqs[j])

            if identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE:
                # union: higher joins lower (lower ID rep)
                if uf.union_lower_wins(i, j):
                    merged += 1

    # Finalize cluster reps
    reps = [uf.find(i) for i in range(n)]

    # Write assignment table: seq_index, fasta_id, rep_index, rep_fasta_id
    assign_path = OUT_DIR / "cluster_assignments.tsv"
    with open(assign_path, "w", newline="\n", encoding="utf-8") as f:
        f.write("seq_index\tfasta_id\trep_index\trep_fasta_id\n")
        for i in range(n):
            r = reps[i]
            f.write(f"{i}\t{ids[i]}\t{r}\t{ids[r]}\n")

    # Write representative IDs
    unique_reps = sorted(set(reps))
    reps_path = OUT_DIR / "nonredundant_reps.txt"
    with open(reps_path, "w", newline="\n", encoding="utf-8") as f:
        for r in unique_reps:
            f.write(str(r) + "\n")

    reps_names_path = OUT_DIR / "nonredundant_rep_names.txt"
    with open(reps_names_path, "w", newline="\n", encoding="utf-8") as f:
        for r in unique_reps:
            f.write(ids[r] + "\n")

    print("\n===== DONE =====")
    print(f"Total sequences: {n}")
    print(f"Pairs tested (after length filter): {tested}")
    print(f"Successful merges: {merged}")
    print(f"Clusters / reps: {len(unique_reps)}")
    print(f"Wrote: {assign_path}")
    print(f"Wrote: {reps_path}")
    print(f"Wrote: {reps_names_path}")


if __name__ == "__main__":
    main()
