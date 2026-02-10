from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, List, Set
from datetime import datetime

from Bio import SeqIO


# ===================== PATHS =====================
BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

# INPUTS
IN_FASTA = OUT / "all_chains_ss_nonzero.fasta"
MMSEQS_TSV = OUT / "mmseqs_clusters_ss_nonzero.tsv"

# Put your real SS counts file name here:
SS_COUNTS = OUT / "ss_counts_20260203_125637.tsv"  # fasta_id\tbeta_res\thelix_res\tbeta_strands\thelices

# OUTPUT TAG (auto so it won't overwrite)
TAG = f"final_nr_by_ss_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

OUT_REP_LIST = OUT / f"nonredundant_rep_names_{TAG}.txt"
OUT_REP_FASTA = OUT / f"nonredundant_rep_chains_{TAG}.fasta"
OUT_DETAILS = OUT / f"cluster_reps_{TAG}.tsv"
OUT_SUMMARY = OUT / f"summary_{TAG}.txt"


# ===================== LOADERS =====================
def load_ss(path: Path) -> Dict[str, Tuple[int, int, int, int]]:
    """
    Expected header:
      fasta_id    beta_res    helix_res    beta_strands    helices
    """
    ss: Dict[str, Tuple[int, int, int, int]] = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline()  # ignore
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
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


def load_mmseqs_rep_member(path: Path) -> Dict[str, Set[str]]:
    """
    MMseqs createtsv gives: rep<TAB>member
    We load rep->set(members), include rep itself.
    """
    clusters: Dict[str, Set[str]] = defaultdict(set)
    total = 0
    bad = 0

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                bad += 1
                continue
            rep, mem = parts[0], parts[1]
            clusters[rep].add(mem)
            total += 1

    # ensure rep in its own set
    for rep in list(clusters.keys()):
        clusters[rep].add(rep)

    return clusters


def build_owner_mapping(raw_clusters: Dict[str, Set[str]], valid_ids: Set[str]) -> Dict[str, str]:
    """
    Critical fix: convert rep->members into a SINGLE owner rep per member.
    Otherwise if your TSV is noisy (or if you later add singletons), you might keep too many.
    Rule: if a member appears under multiple reps, assign it to lexicographically smallest rep.
    """
    owner: Dict[str, str] = {}

    for rep in sorted(raw_clusters.keys()):
        # skip rep not in our fasta universe
        if rep not in valid_ids:
            continue
        for mem in raw_clusters[rep]:
            if mem not in valid_ids:
                continue
            if mem not in owner or rep < owner[mem]:
                owner[mem] = rep

    return owner


def owner_to_clusters(owner: Dict[str, str]) -> Dict[str, List[str]]:
    """
    owner[mem] = rep  -> clusters[rep] = sorted list of members
    """
    clusters: Dict[str, List[str]] = defaultdict(list)
    for mem, rep in owner.items():
        clusters[rep].append(mem)
    for rep in list(clusters.keys()):
        clusters[rep].sort()
    return clusters


# ===================== SELECTION =====================
def pick_best(members: List[str], ss: Dict[str, Tuple[int, int, int, int]]) -> Tuple[str, str]:
    """
    Selection rule (per supervisor message):
    1) If any member has (beta_res + helix_res) > 0, choose max (beta_res + helix_res)
    2) Else choose max (beta_strands + helices)
    tie-break: lexicographically smallest fasta_id
    """
    stats = {m: ss.get(m, (0, 0, 0, 0)) for m in members}

    res_scores = {m: stats[m][0] + stats[m][1] for m in members}
    if any(v > 0 for v in res_scores.values()):
        best = min(members)
        best_sc = -1
        for m in members:
            sc = res_scores[m]
            if sc > best_sc or (sc == best_sc and m < best):
                best_sc = sc
                best = m
        return best, f"residue_sum={best_sc}"

    elem_scores = {m: stats[m][2] + stats[m][3] for m in members}
    best = min(members)
    best_sc = -1
    for m in members:
        sc = elem_scores[m]
        if sc > best_sc or (sc == best_sc and m < best):
            best_sc = sc
            best = m
    return best, f"element_sum={best_sc}"


# ===================== MAIN =====================
def main():
    if not IN_FASTA.exists():
        raise FileNotFoundError(IN_FASTA)
    if not MMSEQS_TSV.exists():
        raise FileNotFoundError(MMSEQS_TSV)
    if not SS_COUNTS.exists():
        raise FileNotFoundError(SS_COUNTS)

    # Load FASTA records (universe of valid ids)
    recs = {r.id: r for r in SeqIO.parse(str(IN_FASTA), "fasta")}
    valid_ids = set(recs.keys())

    # Load SS counts
    ss = load_ss(SS_COUNTS)

    # Load mmseqs rep->members
    raw = load_mmseqs_rep_member(MMSEQS_TSV)

    # Add missing FASTA ids as singleton reps (important: mmseqs TSV may omit singletons)
    appeared = set()
    for rep, mems in raw.items():
        if rep in valid_ids:
            appeared.add(rep)
        for m in mems:
            if m in valid_ids:
                appeared.add(m)

    for fid in sorted(valid_ids):
        if fid not in appeared:
            raw[fid] = {fid}

    # Convert to clean clustering (ONE owner rep per member)
    owner = build_owner_mapping(raw, valid_ids)
    clusters = owner_to_clusters(owner)

    # Select reps
    chosen: List[str] = []
    singletons = 0
    nonsingle = 0

    with OUT_DETAILS.open("w", encoding="utf-8", newline="\n") as fdet:
        fdet.write("cluster_rep\tcluster_size\tchosen_fasta_id\treason\tbeta_res\thelix_res\tbeta_strands\thelices\n")

        for rep in sorted(clusters.keys()):
            members = clusters[rep]
            if not members:
                continue

            if len(members) == 1:
                singletons += 1
                pick = members[0]
                reason = "singleton"
            else:
                nonsingle += 1
                pick, reason = pick_best(members, ss)

            chosen.append(pick)
            br, hr, bs, hx = ss.get(pick, (0, 0, 0, 0))
            fdet.write(f"{rep}\t{len(members)}\t{pick}\t{reason}\t{br}\t{hr}\t{bs}\t{hx}\n")

    # Write outputs
    OUT_REP_LIST.write_text("\n".join(chosen) + "\n", encoding="utf-8")
    SeqIO.write([recs[fid] for fid in chosen], str(OUT_REP_FASTA), "fasta")

    OUT_SUMMARY.write_text(
        "\n".join([
            f"TAG: {TAG}",
            f"Input FASTA chains (ss_nonzero): {len(recs)}",
            f"MMseqs TSV: {MMSEQS_TSV.name}",
            f"SS counts TSV: {SS_COUNTS.name}",
            f"Total clusters (after fixing singletons): {len(clusters)}",
            f"Singleton clusters kept: {singletons}",
            f"Non-singleton clusters: {nonsingle}",
            f"Final representatives: {len(chosen)}",
            "",
            "Wrote:",
            f" - {OUT_REP_LIST.name}",
            f" - {OUT_REP_FASTA.name}",
            f" - {OUT_DETAILS.name}",
            f" - {OUT_SUMMARY.name}",
        ]) + "\n",
        encoding="utf-8"
    )

    print("Wrote:", OUT_REP_LIST)
    print("Wrote:", OUT_REP_FASTA)
    print("Wrote:", OUT_DETAILS)
    print("Wrote:", OUT_SUMMARY)


if __name__ == "__main__":
    main()

# .\venv\Scripts\python.exe .\step3_select_final_reps_by_ss.py
