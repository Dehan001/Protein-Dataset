from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, List, Set, Optional
from Bio import SeqIO


# ===================== PATHS =====================
BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"

# Use your earlier clustering-before-annotation TSV (rep<TAB>member)
MMSEQS_TSV = OUT / "mmseqs_clusters.tsv"

# FASTA that corresponds to that clustering run
IN_FASTA = OUT / "all_chains.fasta"

# Your new SS counts file (update filename!)
SS_COUNTS_TSV = OUT / "ss_counts_20260203_125637.tsv"

# output tag (change to something meaningful)
TAG = "clusterBeforeAnnot_selectBySS"

OUT_REP_LIST = OUT / f"nonredundant_reps_by_ss_{TAG}.txt"
OUT_FASTA = OUT / f"nonredundant_reps_by_ss_{TAG}.fasta"
OUT_SUMMARY = OUT / f"nonredundant_reps_by_ss_{TAG}_summary.txt"
OUT_DETAILS = OUT / f"nonredundant_reps_by_ss_{TAG}_cluster_reps.tsv"
# ===============================================


def load_mmseqs_clusters(path: Path) -> Dict[str, Set[str]]:
    clusters: Dict[str, Set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, mem = parts[0], parts[1]
            clusters[rep].add(mem)
    # ensure rep included
    for rep in list(clusters.keys()):
        clusters[rep].add(rep)
    return clusters


def load_ss_counts(path: Path) -> Dict[str, Tuple[int, int, int, int]]:
    """
    ss_counts format:
      fasta_id  beta_res  helix_res  beta_strands  helices
    """
    ss: Dict[str, Tuple[int, int, int, int]] = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            fid = parts[0].strip()
            try:
                br = int(parts[1])
                hr = int(parts[2])
                bs = int(parts[3])
                hx = int(parts[4])
            except ValueError:
                br = hr = bs = hx = 0
            ss[fid] = (br, hr, bs, hx)
    return ss


def pick_best_member(members: List[str], ss: Dict[str, Tuple[int, int, int, int]]) -> Tuple[str, str]:
    """
    Returns chosen_id, reason
    Priority:
      1) maximize beta_res+helix_res if available (nonzero across the cluster)
      2) else maximize beta_strands+helices
      Tie: lexicographically smallest fasta_id
    """
    # gather stats
    stats = {m: ss.get(m, (0, 0, 0, 0)) for m in members}

    # Determine whether residue counts are "usable" in this cluster
    # (if all members have beta_res+helix_res == 0, we fallback)
    res_scores = {m: (stats[m][0] + stats[m][1]) for m in members}
    use_res = any(v > 0 for v in res_scores.values())

    if use_res:
        # maximize residue score, then tie by id
        best = min(members)  # for tie-breaking
        best_score = -1
        for m in members:
            sc = res_scores[m]
            if sc > best_score or (sc == best_score and m < best):
                best_score = sc
                best = m
        return best, f"residue(beta_res+helix_res)={best_score}"

    # Fallback: maximize counts of elements
    elem_scores = {m: (stats[m][2] + stats[m][3]) for m in members}
    best = min(members)
    best_score = -1
    for m in members:
        sc = elem_scores[m]
        if sc > best_score or (sc == best_score and m < best):
            best_score = sc
            best = m
    return best, f"elements(beta_strands+helices)={best_score}"


def main():
    if not MMSEQS_TSV.exists():
        raise FileNotFoundError(MMSEQS_TSV)
    if not IN_FASTA.exists():
        raise FileNotFoundError(IN_FASTA)
    if not SS_COUNTS_TSV.exists():
        raise FileNotFoundError(SS_COUNTS_TSV)

    # load fasta records
    recs = {r.id: r for r in SeqIO.parse(str(IN_FASTA), "fasta")}

    # clusters and ss
    clusters = load_mmseqs_clusters(MMSEQS_TSV)
    ss = load_ss_counts(SS_COUNTS_TSV)

    # some chains might not appear in MMseqs TSV (depending on createtsv settings)
    all_members = set()
    for rep, mems in clusters.items():
        all_members |= mems

    # treat missing-from-TSV but present-in-FASTA as singleton
    for fid in recs.keys():
        if fid not in all_members and fid not in clusters:
            clusters[fid] = {fid}

    # choose reps
    chosen: List[str] = []
    singletons = 0
    nonsingle = 0

    # details table
    detail_rows: List[Tuple[str, int, str, str, str]] = []
    # rep, size, chosen, reason, top5 scores preview

    for rep in sorted(clusters.keys()):
        members = sorted([m for m in clusters[rep] if m in recs])
        if not members:
            continue

        if len(members) == 1:
            singletons += 1
            chosen_id = members[0]
            chosen.append(chosen_id)
            br, hr, bs, hx = ss.get(chosen_id, (0, 0, 0, 0))
            preview = f"{chosen_id} br={br} hr={hr} bs={bs} hx={hx}"
            detail_rows.append((rep, 1, chosen_id, "singleton", preview))
            continue

        nonsingle += 1
        chosen_id, reason = pick_best_member(members, ss)
        chosen.append(chosen_id)

        # build preview of top few candidates by residue score (then element score)
        scored = []
        for m in members:
            br, hr, bs, hx = ss.get(m, (0, 0, 0, 0))
            scored.append((m, br + hr, bs + hx, br, hr, bs, hx))
        scored.sort(key=lambda x: (-x[1], -x[2], x[0]))
        preview = "; ".join([f"{m} res={r} elem={e} (br={br},hr={hr},bs={bs},hx={hx})"
                             for (m, r, e, br, hr, bs, hx) in scored[:5]])
        detail_rows.append((rep, len(members), chosen_id, reason, preview))

    # write outputs
    OUT_REP_LIST.write_text("\n".join(chosen) + "\n", encoding="utf-8")
    SeqIO.write([recs[fid] for fid in chosen if fid in recs], str(OUT_FASTA), "fasta")

    with OUT_DETAILS.open("w", encoding="utf-8", newline="\n") as f:
        f.write("cluster_rep\tcluster_size\tchosen_fasta_id\treason\ttop_candidates_preview\n")
        for rep, size, chosen_id, reason, preview in detail_rows:
            f.write(f"{rep}\t{size}\t{chosen_id}\t{reason}\t{preview}\n")

    OUT_SUMMARY.write_text(
        f"Input FASTA chains: {len(recs)}\n"
        f"Clusters total: {len(clusters)}\n"
        f"Singleton clusters kept: {singletons}\n"
        f"Non-singleton clusters: {nonsingle}\n"
        f"Final representatives written: {len(chosen)}\n"
        f"SS counts used from: {SS_COUNTS_TSV.name}\n"
        f"MMseqs clusters used from: {MMSEQS_TSV.name}\n"
        f"Outputs:\n"
        f" - {OUT_REP_LIST.name}\n"
        f" - {OUT_FASTA.name}\n"
        f" - {OUT_DETAILS.name}\n",
        encoding="utf-8"
    )

    print("Wrote:", OUT_REP_LIST)
    print("Wrote:", OUT_FASTA)
    print("Wrote:", OUT_DETAILS)
    print("Wrote:", OUT_SUMMARY)


if __name__ == "__main__":
    main()
# .\venv\Scripts\python.exe .\stepE_select_reps_by_ss.py
