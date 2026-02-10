from __future__ import annotations

from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime
from typing import Dict, Set, List, Tuple

# ===================== PATHS =====================
BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

# ---- CHANGE THESE TWO TO YOUR CURRENT RUN ----
FASTA = OUT / "all_chains_ss_nonzero.fasta"              # the FASTA you clustered
MMSEQS_TSV = OUT / "mmseqs_clusters_ss_nonzero.tsv"      # rep<TAB>member from createtsv
# ==============================================

TAG = datetime.now().strftime("%Y%m%d_%H%M%S")

# ===================== OUTPUTS (timestamped) =====================
CHAIN_IDS_OUT = OUT / f"records_chain_ids_{TAG}.txt"

ASSIGN_OUT = OUT / f"records_assignment_{TAG}.tsv"           # fasta_id -> cluster_rep
CLUSTERS_OUT = OUT / f"records_clusters_{TAG}.tsv"           # rep -> size -> members (csv)
REPS_OUT = OUT / f"records_reps_{TAG}.txt"                   # cluster reps (one per cluster)

SINGLETON_CLUSTERS_OUT = OUT / f"records_singleton_clusters_{TAG}.tsv"  # rep -> member
NONSINGLE_CLUSTERS_OUT = OUT / f"records_nonsingle_clusters_{TAG}.tsv"  # rep -> size -> members
SINGLETON_CHAINS_OUT = OUT / f"records_singleton_chains_{TAG}.txt"      # all chains in size=1 clusters
MATCHED_CHAINS_OUT = OUT / f"records_chains_in_size_ge2_{TAG}.txt"      # all chains in size>=2 clusters
MISSING_FROM_TSV_OUT = OUT / f"records_missing_from_tsv_{TAG}.txt"      # chains absent in TSV (added as singleton)

SIZEDIST_OUT = OUT / f"records_cluster_size_distribution_{TAG}.tsv"
SUMMARY_OUT = OUT / f"records_summary_{TAG}.txt"


# ===================== HELPERS =====================
def read_fasta_ids(fasta_path: Path) -> List[str]:
    ids: List[str] = []
    with fasta_path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def load_mmseqs_rep_member(path: Path) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            pairs.append((parts[0], parts[1]))
    return pairs


def main():
    if not FASTA.exists():
        raise FileNotFoundError(f"Missing FASTA: {FASTA}")
    if not MMSEQS_TSV.exists():
        raise FileNotFoundError(f"Missing MMseqs TSV: {MMSEQS_TSV}")

    # ---- Universe from FASTA ----
    all_ids = read_fasta_ids(FASTA)
    id_set = set(all_ids)
    n_total = len(all_ids)

    CHAIN_IDS_OUT.write_text("\n".join(all_ids) + "\n", encoding="utf-8")

    # ---- Read MMseqs pairs ----
    clusters_raw: Dict[str, Set[str]] = defaultdict(set)
    pairs = load_mmseqs_rep_member(MMSEQS_TSV)

    lines_read = 0
    skipped_not_in_fasta = 0
    self_pairs = 0

    for rep, mem in pairs:
        lines_read += 1
        if rep not in id_set or mem not in id_set:
            skipped_not_in_fasta += 1
            continue
        clusters_raw[rep].add(mem)
        if rep == mem:
            self_pairs += 1

    # Ensure rep contains itself
    for rep in list(clusters_raw.keys()):
        clusters_raw[rep].add(rep)

    # ---- Identify who appeared in TSV ----
    appeared: Set[str] = set()
    for rep, mems in clusters_raw.items():
        appeared.add(rep)
        appeared.update(mems)

    # Chains missing from TSV are true singletons (or TSV-generation omitted them)
    missing_from_tsv = [fid for fid in all_ids if fid not in appeared]
    missing_from_tsv_sorted = sorted(missing_from_tsv)
    MISSING_FROM_TSV_OUT.write_text("\n".join(missing_from_tsv_sorted) + "\n", encoding="utf-8")

    # Add missing as singleton clusters
    for fid in missing_from_tsv_sorted:
        clusters_raw[fid] = {fid}

    # ---- Fix ownership in case of weird TSV (member under multiple reps) ----
    # Owner rep chosen as lexicographically smallest rep (deterministic)
    owner_rep: Dict[str, str] = {}
    for rep in sorted(clusters_raw.keys()):
        for mem in clusters_raw[rep]:
            if mem not in owner_rep or rep < owner_rep[mem]:
                owner_rep[mem] = rep

    # Rebuild clean clusters from owner mapping
    clusters: Dict[str, List[str]] = defaultdict(list)
    for mem, rep in owner_rep.items():
        clusters[rep].append(mem)

    for rep in clusters:
        clusters[rep] = sorted(set(clusters[rep]))

    reps_sorted = sorted(clusters.keys())

    # ---- Write assignment ----
    with ASSIGN_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("fasta_id\tcluster_rep\n")
        for fid in all_ids:
            out.write(f"{fid}\t{owner_rep.get(fid, fid)}\n")

    # ---- Write clusters (full) ----
    with CLUSTERS_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("cluster_rep\tcluster_size\tmembers\n")
        for rep in reps_sorted:
            mems = clusters[rep]
            out.write(f"{rep}\t{len(mems)}\t{','.join(mems)}\n")

    # ---- Reps list ----
    REPS_OUT.write_text("\n".join(reps_sorted) + "\n", encoding="utf-8")

    # ---- Singletons vs non-singletons ----
    singleton_chains: List[str] = []
    matched_chains: List[str] = []

    singleton_clusters: List[Tuple[str, str]] = []  # (rep, only_member)
    nonsingle_clusters: List[Tuple[str, int, List[str]]] = []  # (rep, size, members)

    for rep in reps_sorted:
        mems = clusters[rep]
        if len(mems) == 1:
            singleton_chains.append(mems[0])
            singleton_clusters.append((rep, mems[0]))
        else:
            matched_chains.extend(mems)
            nonsingle_clusters.append((rep, len(mems), mems))

    singleton_chains = sorted(set(singleton_chains))
    matched_chains = sorted(set(matched_chains))

    SINGLETON_CHAINS_OUT.write_text("\n".join(singleton_chains) + "\n", encoding="utf-8")
    MATCHED_CHAINS_OUT.write_text("\n".join(matched_chains) + "\n", encoding="utf-8")

    with SINGLETON_CLUSTERS_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("cluster_rep\tonly_member\n")
        for rep, only in singleton_clusters:
            out.write(f"{rep}\t{only}\n")

    with NONSINGLE_CLUSTERS_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("cluster_rep\tcluster_size\tmembers\n")
        for rep, sz, mems in nonsingle_clusters:
            out.write(f"{rep}\t{sz}\t{','.join(mems)}\n")

    # ---- Size distribution ----
    size_counts = Counter(len(clusters[rep]) for rep in reps_sorted)
    with SIZEDIST_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("cluster_size\tnum_clusters\n")
        for s in sorted(size_counts):
            out.write(f"{s}\t{size_counts[s]}\n")

    # ---- Summary ----
    n_clusters = len(reps_sorted)
    n_singleton_clusters = len(singleton_clusters)
    n_nonsingle_clusters = len(nonsingle_clusters)

    # Sanity: chains covered after adding missing should be total
    n_covered = len(owner_rep)

    summary = []
    summary.append("===== MMseqs2 CLUSTER RECORDS SUMMARY =====")
    summary.append(f"FASTA file: {FASTA.name}")
    summary.append(f"MMseqs TSV: {MMSEQS_TSV.name}")
    summary.append("")
    summary.append(f"Total chains in FASTA: {n_total}")
    summary.append(f"TSV lines read: {lines_read}")
    summary.append(f"TSV skipped (IDs not in FASTA): {skipped_not_in_fasta}")
    summary.append(f"TSV self-pairs rep==mem: {self_pairs}")
    summary.append("")
    summary.append(f"Chains missing from TSV (added as singleton clusters): {len(missing_from_tsv_sorted)}")
    summary.append(f"Chains covered by clustering (after adding missing): {n_covered} (should equal {n_total})")
    summary.append("")
    summary.append(f"Total clusters (unique reps): {n_clusters}")
    summary.append(f"Singleton clusters: {n_singleton_clusters}")
    summary.append(f"Non-singleton clusters: {n_nonsingle_clusters}")
    summary.append(f"Chains in singleton clusters: {len(singleton_chains)}")
    summary.append(f"Chains in size>=2 clusters: {len(matched_chains)}")
    summary.append("")
    summary.append("Wrote files:")
    summary.append(f" - {CHAIN_IDS_OUT.name}")
    summary.append(f" - {ASSIGN_OUT.name}")
    summary.append(f" - {CLUSTERS_OUT.name}")
    summary.append(f" - {REPS_OUT.name}")
    summary.append(f" - {SINGLETON_CLUSTERS_OUT.name}")
    summary.append(f" - {NONSINGLE_CLUSTERS_OUT.name}")
    summary.append(f" - {SINGLETON_CHAINS_OUT.name}")
    summary.append(f" - {MATCHED_CHAINS_OUT.name}")
    summary.append(f" - {MISSING_FROM_TSV_OUT.name}")
    summary.append(f" - {SIZEDIST_OUT.name}")

    SUMMARY_OUT.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()

# Wrote: C:\Users\fd02629\Desktop\Protein Experiment\outputs\nonredundant_rep_names_final_nr_by_ss_20260203_223951.txt
# Wrote: C:\Users\fd02629\Desktop\Protein Experiment\outputs\nonredundant_rep_chains_final_nr_by_ss_20260203_223951.fasta
# Wrote: C:\Users\fd02629\Desktop\Protein Experiment\outputs\cluster_reps_final_nr_by_ss_20260203_223951.tsv
# Wrote: C:\Users\fd02629\Desktop\Protein Experiment\outputs\summary_final_nr_by_ss_20260203_223951.txt
# .\venv\Scripts\python.exe .\step_records_mmseqs_clusters.py
