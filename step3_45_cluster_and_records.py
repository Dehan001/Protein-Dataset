from __future__ import annotations
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

FASTA = OUT / "all_chains.fasta"
PAIRS = OUT / "similar_pairs.tsv"

# outputs Dr. Sazzed wants + clustering outputs
CHAIN_IDS_OUT = OUT / "all_chain_ids.txt"
MATCHED_IDS_OUT = OUT / "chains_with_at_least_one_match.txt"
NOMATCH_IDS_OUT = OUT / "chains_with_no_matches.txt"

ASSIGN_OUT = OUT / "cluster_assignments.tsv"   # fasta_id -> rep_fasta_id
CLUSTERS_OUT = OUT / "clusters.tsv"            # rep_fasta_id -> members (comma separated)
REPS_OUT = OUT / "nonredundant_rep_names.txt"  # representative ids (one per cluster)

class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union_lower_wins(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if ra < rb:
            self.parent[rb] = ra
        else:
            self.parent[ra] = rb
        return True

def read_fasta_ids(fasta_path: Path) -> list[str]:
    ids = []
    with fasta_path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def main():
    if not FASTA.exists():
        raise FileNotFoundError(f"Missing: {FASTA}")
    if not PAIRS.exists():
        raise FileNotFoundError(f"Missing: {PAIRS}")

    # ---- chains list (records package) ----
    ids = read_fasta_ids(FASTA)
    n = len(ids)

    CHAIN_IDS_OUT.write_text("\n".join(ids) + "\n", encoding="utf-8")

    # ---- read edges + matched/no-match (records package) ----
    matched = set()
    edges = []
    with PAIRS.open("r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            i = int(parts[0]); j = int(parts[1])
            if i == j:
                continue
            edges.append((i, j))
            matched.add(i); matched.add(j)

    matched_list = sorted(matched)
    nomatch_list = [i for i in range(n) if i not in matched]

    MATCHED_IDS_OUT.write_text("\n".join(ids[i] for i in matched_list) + "\n", encoding="utf-8")
    NOMATCH_IDS_OUT.write_text("\n".join(ids[i] for i in nomatch_list) + "\n", encoding="utf-8")

    # ---- Step 4 clustering (lower-id wins) ----
    uf = UnionFind(n)
    merges = 0
    for (i, j) in edges:
        if uf.union_lower_wins(i, j):
            merges += 1

    reps = [uf.find(i) for i in range(n)]

    # group members by representative
    clusters = defaultdict(list)
    for i in range(n):
        clusters[reps[i]].append(i)

    rep_indices_sorted = sorted(clusters.keys())

    # write cluster_assignments.tsv using FASTA IDs
    with ASSIGN_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("seq_index\tfasta_id\trep_index\trep_fasta_id\n")
        for i in range(n):
            r = reps[i]
            out.write(f"{i}\t{ids[i]}\t{r}\t{ids[r]}\n")

    # write clusters.tsv (rep id -> member ids)
    with CLUSTERS_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("rep_fasta_id\tmembers\n")
        for r in rep_indices_sorted:
            member_ids = ",".join(ids[m] for m in clusters[r])
            out.write(f"{ids[r]}\t{member_ids}\n")

    # representatives list (placeholder selection = rep id = lowest id in cluster)
    REPS_OUT.write_text("\n".join(ids[r] for r in rep_indices_sorted) + "\n", encoding="utf-8")

    # ---- print summary for email/paper ----
    print("===== RECORDS + CLUSTERING SUMMARY =====")
    print(f"Total chains (FASTA records): {n}")
    print(f"Pairs/edges passing thresholds (rows in similar_pairs.tsv): {len(edges)}")
    print(f"Chains with >=1 match: {len(matched_list)}")
    print(f"Chains with no matches: {len(nomatch_list)}")
    print(f"Clusters (unique representatives): {len(rep_indices_sorted)}")
    print(f"Union merges performed: {merges}")
    print("Wrote:")
    print(" -", CHAIN_IDS_OUT)
    print(" -", MATCHED_IDS_OUT)
    print(" -", NOMATCH_IDS_OUT)
    print(" -", ASSIGN_OUT)
    print(" -", CLUSTERS_OUT)
    print(" -", REPS_OUT)

if __name__ == "__main__":
    main()
