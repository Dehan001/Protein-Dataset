from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, List
from collections import defaultdict
from Bio import SeqIO

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "outputs"

FASTA_PATH = OUT_DIR / "all_chains.fasta"
BETA_COUNTS_TSV = OUT_DIR / "beta_counts.tsv"

# Choose ONE input cluster file:
MMSEQS_TSV = OUT_DIR / "mmseqs_clusters.tsv"   # rep_id <tab> member_id (one per line)
CLUSTERS_TSV = OUT_DIR / "clusters.tsv"        # rep_index <tab> members (comma-separated indices)  [optional]

REP_OUT = OUT_DIR / "nonredundant_rep_names_by_beta.txt"
# ====================================================


def load_beta_counts(path: Path) -> Dict[Tuple[str, str], int]:
    beta: Dict[Tuple[str, str], int] = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            pdb_id = parts[0].strip().upper()
            chain = parts[1].strip()
            try:
                cnt = int(parts[2])
            except ValueError:
                continue
            beta[(pdb_id, chain)] = cnt
    return beta


def parse_fasta_id(fid: str) -> Tuple[str, str]:
    """
    Expect: PDBID_EMD-XXXXX_CHAIN
    We extract pdb_id and chain_id (last token).
    """
    parts = fid.split("_")
    pdb_id = parts[0].upper()
    chain = parts[-1]
    return pdb_id, chain


def load_clusters_from_mmseqs(tsv_path: Path) -> Dict[str, List[str]]:
    """
    mmseqs createtsv produces lines: rep_id \t member_id
    We build rep -> [members...]
    """
    clusters: Dict[str, List[str]] = defaultdict(list)
    with tsv_path.open("r", encoding="utf-8") as f:
        for line in f:
            rep, mem = line.rstrip("\n").split("\t")[:2]
            clusters[rep].append(mem)
    return clusters


def load_clusters_from_indices(clusters_path: Path, ids: List[str]) -> Dict[str, List[str]]:
    """
    clusters.tsv: rep_index \t members(comma sep indices)
    Convert indices -> fasta IDs.
    """
    clusters: Dict[str, List[str]] = {}
    with clusters_path.open("r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            rep_idx_str, members_str = line.rstrip("\n").split("\t")[:2]
            member_idxs = [int(x) for x in members_str.split(",") if x.strip()]
            rep_id = ids[int(rep_idx_str)]
            clusters[rep_id] = [ids[i] for i in member_idxs]
    return clusters


def main():
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f"Missing {FASTA_PATH}")
    if not BETA_COUNTS_TSV.exists():
        raise FileNotFoundError(f"Missing {BETA_COUNTS_TSV} (run Step 5.1 first)")

    records = list(SeqIO.parse(str(FASTA_PATH), "fasta"))
    ids = [r.id for r in records]
    id_set = set(ids)

    beta = load_beta_counts(BETA_COUNTS_TSV)

    # Load clusters
    clusters: Dict[str, List[str]]
    if MMSEQS_TSV.exists():
        clusters = load_clusters_from_mmseqs(MMSEQS_TSV)
        # sanity: keep only IDs that exist in FASTA
        clusters = {rep: [m for m in mems if m in id_set] for rep, mems in clusters.items() if rep in id_set}
        source = "mmseqs_clusters.tsv"
    elif CLUSTERS_TSV.exists():
        clusters = load_clusters_from_indices(CLUSTERS_TSV, ids)
        source = "clusters.tsv"
    else:
        raise FileNotFoundError("Provide either outputs/mmseqs_clusters.tsv or outputs/clusters.tsv")

    reps: List[str] = []

    # Choose max-beta representative per cluster (tie-break: lexicographically smallest ID)
    for rep, members in clusters.items():
        best_id = None
        best_beta = -1

        for fid in members:
            pdb_id, chain = parse_fasta_id(fid)
            b = beta.get((pdb_id, chain), 0)
            if b > best_beta or (b == best_beta and (best_id is None or fid < best_id)):
                best_beta = b
                best_id = fid

        if best_id is None:
            # fallback: keep rep
            best_id = rep

        reps.append(best_id)

    reps = sorted(set(reps))
    REP_OUT.write_text("\n".join(reps) + "\n", encoding="utf-8")

    print("===== STEP 5.2 DONE =====")
    print("Cluster source:", source)
    print("Clusters:", len(clusters))
    print("Selected reps:", len(reps))
    print("Wrote:", REP_OUT)


if __name__ == "__main__":
    main()

# .\venv\Scripts\python.exe .\step5_2_select_reps_by_beta.py
