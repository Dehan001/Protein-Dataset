from __future__ import annotations

from pathlib import Path
from collections import defaultdict, Counter

# ------------------- PATHS -------------------
BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

FASTA = OUT / "all_chains.fasta"

# Put your MMseqs TSV here (rep \t member). Example file name:
#   outputs/mmseqs_clusters.tsv
# or whatever you created from `mmseqs createtsv ...`
MMSEQS_TSV = OUT / "mmseqs_clusters.tsv"

# ------------------- OUTPUTS (records package) -------------------
CHAIN_IDS_OUT = OUT / "all_chain_ids.txt"

# clustering outputs
ASSIGN_OUT = OUT / "mmseqs_cluster_assignments.tsv"   # fasta_id -> rep_fasta_id
CLUSTERS_OUT = OUT / "mmseqs_clusters.tsv"            # rep_fasta_id -> members (comma sep)
REPS_OUT = OUT / "mmseqs_nonredundant_rep_names.txt"  # list of reps
SINGLETONS_OUT = OUT / "mmseqs_singletons.txt"        # chains with no matches (singleton clusters)
MATCHED_OUT = OUT / "mmseqs_chains_with_at_least_one_match.txt"  # in cluster size>=2

# summary stats
SUMMARY_OUT = OUT / "mmseqs_records_summary.txt"
SIZEDIST_OUT = OUT / "mmseqs_cluster_size_distribution.tsv"


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
            f"-> Put your mmseqs TSV at that path, or change MMSEQS_TSV in the script."
        )

    # ---- total chain list (paper record) ----
    all_ids = read_fasta_ids(FASTA)
    n_total = len(all_ids)
    id_set = set(all_ids)

    CHAIN_IDS_OUT.write_text("\n".join(all_ids) + "\n", encoding="utf-8")

    # ---- read MMseqs rep->member pairs ----
    # Build clusters: rep_id -> set(members)
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

            # optional sanity: ignore IDs not in our FASTA universe
            if rep not in id_set or mem not in id_set:
                skipped += 1
                continue

            clusters[rep].add(mem)
            if rep == mem:
                self_pairs += 1

    # Some MMseqs outputs can omit singletons depending on how you generated TSV.
    # If a chain never appears anywhere in the TSV, it is effectively a singleton.
    appeared = set()
    for rep, mems in clusters.items():
        appeared.add(rep)
        appeared.update(mems)

    # Ensure every chain is accounted for:
    # - If a chain never appeared in TSV, it is a singleton cluster by itself.
    for fid in all_ids:
        if fid not in appeared:
            clusters[fid].add(fid)

    # Ensure each rep contains itself (just in case TSV missed self-line)
    for rep in list(clusters.keys()):
        clusters[rep].add(rep)

    # ---- normalize: some members might appear under multiple reps in weird TSVs ----
    # MMseqs output should be a proper clustering, but we’ll enforce a consistent assignment:
    # - If a chain appears in multiple clusters, pick the lexicographically smallest rep as owner.
    #   (You can change to "first seen" if you prefer.)
    owner_rep: dict[str, str] = {}
    for rep in sorted(clusters.keys()):
        for mem in clusters[rep]:
            if mem not in owner_rep:
                owner_rep[mem] = rep
            else:
                # already assigned; keep smallest rep
                if rep < owner_rep[mem]:
                    owner_rep[mem] = rep

    # rebuild clusters from owner mapping
    fixed_clusters: dict[str, list[str]] = defaultdict(list)
    for mem, rep in owner_rep.items():
        fixed_clusters[rep].append(mem)

    # sort members for stable output
    for rep in fixed_clusters:
        fixed_clusters[rep] = sorted(fixed_clusters[rep])

    reps_sorted = sorted(fixed_clusters.keys())

    # ---- write assignment table ----
    with ASSIGN_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("fasta_id\trep_fasta_id\n")
        for fid in all_ids:
            rep = owner_rep.get(fid, fid)
            out.write(f"{fid}\t{rep}\n")

    # ---- write clusters.tsv ----
    with CLUSTERS_OUT.open("w", encoding="utf-8", newline="\n") as out:
        out.write("rep_fasta_id\tcluster_size\tmembers\n")
        for rep in reps_sorted:
            members = fixed_clusters[rep]
            out.write(f"{rep}\t{len(members)}\t{','.join(members)}\n")

    # ---- representatives list ----
    REPS_OUT.write_text("\n".join(reps_sorted) + "\n", encoding="utf-8")

    # ---- compute “match vs no match” (paper record) ----
    # Define:
    # - "no matches" = singleton clusters (size == 1)
    # - ">=1 match"  = clusters size >= 2, and all members inside them
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

    # ---- summary for emailing/paper ----
    n_clusters = len(reps_sorted)
    n_singletons = len(singletons)
    n_matched = len(matched)
    n_covered = len(owner_rep)  # should equal n_total after our "account for missing" step

    summary = []
    summary.append("===== MMseqs2 RECORDS SUMMARY =====")
    summary.append(f"FASTA total chains: {n_total}")
    summary.append(f"MMseqs TSV lines read: {lines}")
    summary.append(f"MMseqs TSV skipped (bad/missing IDs): {skipped}")
    summary.append(f"MMseqs self-pairs (rep==mem) seen: {self_pairs}")
    summary.append(f"Chains covered by clustering (after adding missing as singletons): {n_covered}")
    summary.append(f"Total clusters (representatives): {n_clusters}")
    summary.append(f"Singleton clusters (no matches): {n_singletons}")
    summary.append(f"Chains with >=1 match (in cluster size>=2): {n_matched}")
    summary.append("")
    summary.append("Wrote:")
    summary.append(f" - {CHAIN_IDS_OUT}")
    summary.append(f" - {ASSIGN_OUT}")
    summary.append(f" - {CLUSTERS_OUT}")
    summary.append(f" - {REPS_OUT}")
    summary.append(f" - {SINGLETONS_OUT}")
    summary.append(f" - {MATCHED_OUT}")
    summary.append(f" - {SIZEDIST_OUT}")

    SUMMARY_OUT.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
