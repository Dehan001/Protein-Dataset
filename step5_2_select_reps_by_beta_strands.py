from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

# ============== SETTINGS ==============
BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"

MMSEQS_CLUSTERS = OUT / "mmseqs_clusters.tsv"     # rep<TAB>member
BETA_COUNTS = OUT / "beta_counts.tsv"             # fasta_id<TAB>beta_strands
OUT_REPS = OUT / "nonredundant_rep_names.txt"     # one fasta_id per line
OUT_CLUSTER_REPS = OUT / "cluster_reps.tsv"       # rep_id<TAB>chosen_id<TAB>beta
# =====================================


def load_beta(path: Path) -> Dict[str, int]:
    beta: Dict[str, int] = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            fid = parts[0]
            try:
                b = int(parts[1])
            except ValueError:
                continue
            beta[fid] = b
    return beta


def load_clusters(path: Path) -> Dict[str, List[str]]:
    """
    mmseqs createtsv typically gives rep -> member lines.
    We'll build: rep -> [members...]
    """
    clusters: Dict[str, List[str]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rep, mem = line.split("\t")[:2]
            clusters[rep].append(mem)
    return clusters


def main():
    if not MMSEQS_CLUSTERS.exists():
        raise FileNotFoundError(f"Missing {MMSEQS_CLUSTERS}")
    if not BETA_COUNTS.exists():
        raise FileNotFoundError(f"Missing {BETA_COUNTS}. Run step5_1_make_beta_counts_dssp.py first.")

    beta = load_beta(BETA_COUNTS)
    clusters = load_clusters(MMSEQS_CLUSTERS)

    chosen: List[str] = []
    cluster_lines: List[Tuple[str, str, int]] = []

    # Select: max beta_strands; tie-break: lexicographically smallest fasta_id
    for rep, members in sorted(clusters.items()):
        best_id = None
        best_beta = -1

        for m in members:
            b = beta.get(m, -1)
            if b > best_beta:
                best_beta = b
                best_id = m
            elif b == best_beta and best_id is not None and m < best_id:
                best_id = m

        # If everything missing beta (e.g., DSSP failed), fall back to rep itself
        if best_id is None:
            best_id = rep
            best_beta = beta.get(rep, -1)

        chosen.append(best_id)
        cluster_lines.append((rep, best_id, best_beta))

    OUT_REPS.write_text("\n".join(chosen) + "\n", encoding="utf-8")

    with OUT_CLUSTER_REPS.open("w", encoding="utf-8", newline="\n") as out:
        out.write("mmseqs_rep\tchosen_rep\tbeta_strands\n")
        for rep, best_id, b in cluster_lines:
            out.write(f"{rep}\t{best_id}\t{b}\n")

    print("===== STEP 5.2 DONE =====")
    print("Clusters:", len(clusters))
    print("Wrote:", OUT_REPS)
    print("Wrote:", OUT_CLUSTER_REPS)
    print("Note: If DSSP failed for some chains, their beta_strands may be -1 and selection falls back safely.")


if __name__ == "__main__":
    main()

# .\venv\Scripts\python.exe .\step5_2_select_reps_by_beta_strands.py
