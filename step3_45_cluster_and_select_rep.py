from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


# ---------------- Paths ----------------
BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)

SIMILAR_PAIRS = OUT_DIR / "similar_pairs.tsv"          # from Step 1–2 NW
BETA_COUNTS = OUT_DIR / "beta_counts.tsv"              # optional (Step 5)
ASSIGN_OUT = OUT_DIR / "cluster_assignments.tsv"
CLUSTERS_OUT = OUT_DIR / "clusters.tsv"
REP_OUT = OUT_DIR / "nonredundant_rep_names.txt"
# --------------------------------------


class UnionFind:
    """Union-Find where lower representative always wins (higher attaches to lower)."""
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


def read_beta_counts(path: Path) -> Optional[Dict[str, int]]:
    """
    Expected TSV:
    fasta_id    beta_strands
    9CTH_EMD-42405_B    17
    """
    if not path.exists():
        return None

    beta: Dict[str, int] = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 2:
            raise ValueError("beta_counts.tsv must have at least 2 columns: fasta_id<TAB>beta_strands")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            fid = parts[0].strip()
            try:
                beta[fid] = int(parts[1])
            except ValueError:
                continue
    return beta


def main():
    if not SIMILAR_PAIRS.exists():
        raise FileNotFoundError(f"Missing {SIMILAR_PAIRS}. Run Step 1–2 first.")

    # We will parse BOTH numeric indices and fasta IDs from similar_pairs.tsv
    # Your Step1–2 file header is typically:
    # i j id_i id_j identity coverage_shorter aligned_res len_shorter len_i len_j
    with SIMILAR_PAIRS.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: idx for idx, name in enumerate(header)}

        required = ["i", "j", "id_i", "id_j"]
        missing = [c for c in required if c not in col]
        if missing:
            raise ValueError(
                f"similar_pairs.tsv is missing columns {missing}. "
                f"Found columns: {header}"
            )

        edges: List[Tuple[int, int]] = []
        id_by_index: Dict[int, str] = {}
        max_idx = -1

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(col.values()):
                continue

            i = int(parts[col["i"]])
            j = int(parts[col["j"]])
            if i == j:
                continue

            id_i = parts[col["id_i"]].strip()
            id_j = parts[col["id_j"]].strip()

            id_by_index[i] = id_i
            id_by_index[j] = id_j

            if i > max_idx: max_idx = i
            if j > max_idx: max_idx = j

            edges.append((i, j))

    if max_idx < 0:
        raise ValueError("similar_pairs.tsv has no data rows / edges.")

    n = max_idx + 1
    uf = UnionFind(n)

    merges = 0
    for (i, j) in edges:
        if uf.union_lower_wins(i, j):
            merges += 1

    reps = [uf.find(i) for i in range(n)]

    # Group members by rep_index
    clusters: Dict[int, List[int]] = defaultdict(list)
    for idx in range(n):
        clusters[reps[idx]].append(idx)

    rep_indices = sorted(clusters.keys())

    # ---- Write cluster_assignments.tsv (with fasta IDs) ----
    with ASSIGN_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("seq_index\tfasta_id\trep_index\trep_fasta_id\n")
        for i in range(n):
            r = reps[i]
            out.write(f"{i}\t{id_by_index.get(i,'')}\t{r}\t{id_by_index.get(r,'')}\n")

    # ---- Write clusters.tsv (rep_fasta_id -> members_fasta_ids) ----
    with CLUSTERS_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("rep_index\trep_fasta_id\tmembers\n")
        for r in rep_indices:
            rep_fid = id_by_index.get(r, "")
            members = [id_by_index.get(m, "") for m in clusters[r]]
            members = [m for m in members if m]  # drop blanks
            out.write(f"{r}\t{rep_fid}\t{','.join(members)}\n")

    # ---- Step 5: choose representative per cluster ----
    beta = read_beta_counts(BETA_COUNTS)

    chosen_rep_fasta_ids: List[str] = []

    for r in rep_indices:
        members_idx = clusters[r]
        members_fid = [(m, id_by_index.get(m, "")) for m in members_idx]
        members_fid = [(m, fid) for (m, fid) in members_fid if fid]

        if not members_fid:
            continue

        if beta is None:
            # Placeholder: keep the Step 4 representative (lowest-ID cluster root)
            chosen_rep_fasta_ids.append(id_by_index.get(r, ""))
            continue

        # Choose member with max beta_strands; tie-break by lower index
        best_idx, best_fid = members_fid[0]
        best_beta = beta.get(best_fid, -1)

        for m_idx, m_fid in members_fid[1:]:
            b = beta.get(m_fid, -1)
            if b > best_beta or (b == best_beta and m_idx < best_idx):
                best_idx, best_fid, best_beta = m_idx, m_fid, b

        chosen_rep_fasta_ids.append(best_fid)

    REP_OUT.write_text("\n".join([x for x in chosen_rep_fasta_ids if x]) + "\n", encoding="utf-8")

    print("===== STEP 4–5 DONE (from similar_pairs.tsv) =====")
    print("Total nodes:", n)
    print("Edges used:", len(edges))
    print("Union merges performed:", merges)
    print("Clusters:", len(rep_indices))
    print("Wrote:", ASSIGN_OUT)
    print("Wrote:", CLUSTERS_OUT)
    print("Wrote:", REP_OUT)
    if beta is None:
        print("NOTE: beta_counts.tsv not found -> reps = Step 4 roots (lowest-ID).")
    else:
        print("beta_counts.tsv found -> reps chosen by max beta_strands per cluster.")


if __name__ == "__main__":
    main()

# .\.venv\Scripts\python.exe .\step3_45_cluster_and_select_rep.py
