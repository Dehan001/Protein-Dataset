from __future__ import annotations

from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime

# ------------------- PATHS -------------------
BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

# CURRENT inputs (beta-only pipeline)
FASTA = OUT / "all_chains_beta_only.fasta"
MMSEQS_TSV = OUT / "mmseqs_clusters_beta_only.tsv"   # rep \t member

# ------------------- OUTPUT FOLDER (no overwrite) -------------------
# Creates: outputs/records_mmseqs_beta_only_YYYYMMDD_HHMMSS/
RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S")
REC_DIR = OUT / f"records_mmseqs_beta_only_{RUN_ID}"

def ensure_new_dir(p: Path):
    if p.exists():
        raise FileExistsError(f"Refusing to overwrite existing folder: {p}")
    p.mkdir(parents=True, exist_ok=False)

# ------------------- OUTPUTS -------------------
CHAIN_IDS_OUT   = REC_DIR / "all_chain_ids_beta_only.txt"

ASSIGN_OUT      = REC_DIR / "mmseqs_cluster_assignments_beta_only.tsv"   # fasta_id -> rep_fasta_id
CLUSTERS_OUT    = REC_DIR / "mmseqs_clusters_beta_only_collapsed.tsv"    # rep -> members (comma sep)
REPS_OUT        = REC_DIR / "mmseqs_nonredundant_rep_names_beta_only.txt"
SINGLETONS_OUT  = REC_DIR / "mmseqs_singletons_beta_only.txt"
MATCHED_OUT     = REC_DIR / "mmseqs_chains_with_at_least_one_match_beta_only.txt"

SUMMARY_OUT     = REC_DIR / "mmseqs_records_summary_beta_only.txt"
SIZEDIST_OUT    = REC_DIR / "mmseqs_cluster_size_distribution_beta_only.tsv"


def read_fasta_ids(fasta_path: Path) -> list[str]:
    ids = []
    with fasta_path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def main():
    if not FASTA.exists():
        raise FileNotFoundError(f"Missing FASTA: {FASTA}")
    if not MMSEQS_TSV.exists():
        raise FileNotFoundError(
            f"Missing MMseqs TSV: {MMSEQS_TSV}\n"
            f"-> Make sure you ran createtsv and wrote it to that filename."
        )

    ensure_new_dir(REC_DIR)

    # ---- total chain list (paper record) ----
    all_ids = read_fasta_ids(FASTA)
    n_total = len(all_ids)
    id_set = set(all_ids)

    CHAIN_IDS_OUT.write_text("\n".join(all_ids) + "\n", encoding="utf-8")

    # ---- read MMseqs rep->member pairs ----
    clusters: dict[str, set[str]] = defaultdict(set)

    lines = 0
    skipped = 0
    self_pairs = 0

    with MMSEQS_TSV.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                skipped += 1
                continue
            rep, mem = parts[0], parts[1]
            lines += 1

            # ignore IDs not in our FASTA universe
            if rep not in id_set or mem not in id_set:
                skipped += 1
                continue

            clusters[rep].add(mem)
            if rep == mem:
                self_pairs += 1

    # ---- add missing as singletons (MMseqs TSV may omit them) ----
    appeared = set()
    for rep, mems in clusters.items():
        appeared.add(rep)
        appeared.update(mems)

    for fid in all_ids:
        if fid not in appeared:
            clusters[fid].add(fid)

    for rep in list(clusters.keys()):
        clusters[rep].add(rep)

    # ---- enforce consistent ownership (rare multi-rep overlaps) ----
    owner_rep: dict[str, str] = {}
    for rep in sorted(clusters.keys()):
        for mem in clusters[rep]:
            if mem not in owner_rep:
                owner_rep[mem] = rep
            else:
                if rep < owner_rep[mem]:
                    owner_rep[mem] = rep

    fixed_clusters: dict[str, list[str]] = defaultdict(list)
    for mem, rep in owner_rep.items():
        fixed_clusters[rep].append(mem)

    for rep in fixed_clusters:
        fixed_clusters[rep] = sorted(fixed_clusters[rep])

    reps_sorted = sorted(fixed_clusters.keys())

    # ---- write assignment table ----
    with ASSIGN_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("fasta_id\trep_fasta_id\n")
        for fid in all_ids:
            rep = owner_rep.get(fid, fid)
            out.write(f"{fid}\t{rep}\n")

    # ---- write collapsed clusters ----
    with CLUSTERS_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("rep_fasta_id\tcluster_size\tmembers\n")
        for rep in reps_sorted:
            members = fixed_clusters[rep]
            out.write(f"{rep}\t{len(members)}\t{','.join(members)}\n")

    # ---- representatives list ----
    REPS_OUT.write_text("\n".join(reps_sorted) + "\n", encoding="utf-8")

    # ---- singleton vs matched ----
    singletons = []
    matched = []
    for rep in reps_sorted:
        members = fixed_clusters[rep]
        if len(members) == 1:
            singletons.append(members[0])
        else:
            matched.extend(members)

    singletons = sorted(set(singletons))
    matched = sorted(set(matched))

    SINGLETONS_OUT.write_text("\n".join(singletons) + "\n", encoding="utf-8")
    MATCHED_OUT.write_text("\n".join(matched) + "\n", encoding="utf-8")

    # ---- cluster size distribution ----
    size_counts = Counter(len(fixed_clusters[rep]) for rep in reps_sorted)
    with SIZEDIST_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("cluster_size\tnum_clusters\n")
        for size in sorted(size_counts):
            out.write(f"{size}\t{size_counts[size]}\n")

    # ---- summary ----
    n_clusters = len(reps_sorted)
    n_singletons = len(singletons)
    n_matched = len(matched)
    n_covered = len(owner_rep)

    summary = []
    summary.append("===== MMseqs2 RECORDS SUMMARY (BETA-ONLY) =====")
    summary.append(f"Run folder: {REC_DIR}")
    summary.append(f"Input FASTA: {FASTA}")
    summary.append(f"Input MMseqs TSV: {MMSEQS_TSV}")
    summary.append("")
    summary.append(f"FASTA total chains: {n_total}")
    summary.append(f"MMseqs TSV lines read: {lines}")
    summary.append(f"MMseqs TSV skipped (bad/missing IDs): {skipped}")
    summary.append(f"MMseqs self-pairs (rep==mem) seen: {self_pairs}")
    summary.append(f"Chains covered (after adding missing as singletons): {n_covered}")
    summary.append(f"Total clusters (representatives): {n_clusters}")
    summary.append(f"Singleton clusters: {n_singletons}")
    summary.append(f"Chains in clusters size>=2: {n_matched}")
    summary.append("")
    summary.append("Wrote:")
    summary.append(f" - {CHAIN_IDS_OUT.name}")
    summary.append(f" - {ASSIGN_OUT.name}")
    summary.append(f" - {CLUSTERS_OUT.name}")
    summary.append(f" - {REPS_OUT.name}")
    summary.append(f" - {SINGLETONS_OUT.name}")
    summary.append(f" - {MATCHED_OUT.name}")
    summary.append(f" - {SIZEDIST_OUT.name}")
    summary.append(f" - {SUMMARY_OUT.name}")

    SUMMARY_OUT.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
# .\venv\Scripts\python.exe .\step_records_mmseqs_beta_only.py
