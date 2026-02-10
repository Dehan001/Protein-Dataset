from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, Set, List, Tuple
from Bio import SeqIO

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"

FASTA = OUT / "all_chains_beta_only.fasta"
BETA_TSV = OUT / "beta_counts.tsv"  # 2 columns: fasta_id, beta_strands
MMSEQS_TSV = OUT / "mmseqs_clusters_beta_only.tsv"  # rep<TAB>member from createtsv

# New output folder (no overwrite)
RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S")
DEST = OUT / f"final_reps_beta_only_{RUN_ID}"
DEST.mkdir(parents=True, exist_ok=False)

OUT_REP_LIST = DEST / "nonredundant_rep_names.txt"
OUT_FASTA = DEST / "nonredundant_rep_chains.fasta"
OUT_CLUSTER_REPS = DEST / "cluster_reps.tsv"
OUT_SUMMARY = DEST / "summary.txt"


def load_beta_counts_2col(path: Path) -> Dict[str, int]:
    beta: Dict[str, int] = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header != ["fasta_id", "beta_strands"]:
            raise ValueError("beta_counts.tsv must be: fasta_id<TAB>beta_strands (2 columns)")
        for line in f:
            if not line.strip():
                continue
            fid, b = line.rstrip("\n").split("\t")
            try:
                beta[fid] = int(b)
            except ValueError:
                beta[fid] = 0
    return beta


def load_mmseqs_clusters(path: Path) -> Dict[str, Set[str]]:
    # createtsv gives rep\tmember lines (often includes rep->rep)
    clusters: Dict[str, Set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rep, mem = line.split("\t")[:2]
            clusters[rep].add(mem)
    # ensure rep included
    for rep in list(clusters.keys()):
        clusters[rep].add(rep)
    return clusters


def main():
    if not FASTA.exists():
        raise FileNotFoundError(FASTA)
    if not BETA_TSV.exists():
        raise FileNotFoundError(BETA_TSV)
    if not MMSEQS_TSV.exists():
        raise FileNotFoundError(MMSEQS_TSV)

    recs = {r.id: r for r in SeqIO.parse(str(FASTA), "fasta")}
    beta = load_beta_counts_2col(BETA_TSV)
    clusters = load_mmseqs_clusters(MMSEQS_TSV)

    # Add any fasta ids missing from mmseqs TSV as singleton clusters
    all_members = set()
    for rep, mems in clusters.items():
        all_members |= mems

    for fid in recs.keys():
        if fid not in all_members and fid not in clusters:
            clusters[fid] = {fid}

    # Select reps
    chosen_final: List[str] = []
    singleton_clusters = 0
    nonsingle_clusters = 0

    with OUT_CLUSTER_REPS.open("w", encoding="utf-8", newline="\n") as out:
        out.write("cluster_rep\tchosen_fasta_id\tbeta_strands\tcluster_size\n")

        for rep in sorted(clusters.keys()):
            members = sorted([m for m in clusters[rep] if m in recs])
            if not members:
                continue

            if len(members) == 1:
                singleton_clusters += 1
                chosen = members[0]
                chosen_final.append(chosen)
                out.write(f"{rep}\t{chosen}\t{beta.get(chosen,0)}\t1\n")
                continue

            nonsingle_clusters += 1

            # Choose highest beta_strands; tie-break by lexicographically smallest fasta_id (deterministic)
            best = None
            best_b = -1
            for m in members:
                b = beta.get(m, 0)
                if b > best_b:
                    best_b = b
                    best = m
                elif b == best_b and best is not None and m < best:
                    best = m

            chosen_final.append(best)
            out.write(f"{rep}\t{best}\t{best_b}\t{len(members)}\n")

    # Write final list + fasta
    OUT_REP_LIST.write_text("\n".join(chosen_final) + "\n", encoding="utf-8")
    SeqIO.write([recs[fid] for fid in chosen_final if fid in recs], str(OUT_FASTA), "fasta")

    OUT_SUMMARY.write_text(
        f"Input FASTA chains (beta-only): {len(recs)}\n"
        f"Total clusters (after adding missing as singletons): {len(clusters)}\n"
        f"Singleton clusters kept: {singleton_clusters}\n"
        f"Non-singleton clusters: {nonsingle_clusters}\n"
        f"Final representatives written: {len(chosen_final)}\n"
        f"Output folder: {DEST}\n",
        encoding="utf-8"
    )

    print("DONE")
    print("Output folder:", DEST)
    print("Singleton clusters:", singleton_clusters)
    print("Non-singleton clusters:", nonsingle_clusters)
    print("Final reps:", len(chosen_final))


if __name__ == "__main__":
    main()
