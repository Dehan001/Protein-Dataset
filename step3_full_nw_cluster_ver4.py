from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
import multiprocessing as mp
import time

# -------- SETTINGS --------
BASE = Path(__file__).resolve().parent
FASTA_PATH = BASE / "outputs" / "all_chains.fasta"
MMSEQS_TSV = BASE / "outputs" / "mmseqs_clusters.tsv"
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)

MIN_IDENTITY = 0.90
MIN_COVERAGE = 0.80

GAP_OPEN = -10.0
GAP_EXTEND = -0.5

# Use 8-10 as requested
NW_WORKERS = 8  # change to 10 if you want
PAIRS_PER_CHUNK = 3000
# -------------------------

# Worker globals
G_SEQS: List[str] | None = None
G_IDS: List[str] | None = None
G_LENS: List[int] | None = None
G_BLOSUM62 = None
G_GAP_OPEN: float | None = None
G_GAP_EXTEND: float | None = None

def _init_worker(seqs: List[str], ids: List[str], lens: List[int], gap_open: float, gap_extend: float):
    global G_SEQS, G_IDS, G_LENS, G_BLOSUM62, G_GAP_OPEN, G_GAP_EXTEND
    G_SEQS = seqs
    G_IDS = ids
    G_LENS = lens
    G_GAP_OPEN = gap_open
    G_GAP_EXTEND = gap_extend
    G_BLOSUM62 = substitution_matrices.load("BLOSUM62")

def _nw_stats(seq1: str, seq2: str) -> Tuple[float, float, int, int]:
    """
    identity = matches / aligned_residues
    coverage_shorter = aligned_residues / len(shorter chain)
    aligned_residues = positions where both are not gaps
    """
    alns = pairwise2.align.globalds(seq1, seq2, G_BLOSUM62, G_GAP_OPEN, G_GAP_EXTEND, one_alignment_only=True)
    if not alns:
        shorter = min(len(seq1), len(seq2))
        return 0.0, 0.0, 0, shorter

    a1, a2, score, start, end = alns[0]
    matches = 0
    aligned = 0
    for c1, c2 in zip(a1, a2):
        if c1 != "-" and c2 != "-":
            aligned += 1
            if c1 == c2:
                matches += 1

    shorter = min(len(seq1), len(seq2))
    if aligned == 0 or shorter == 0:
        return 0.0, 0.0, 0, shorter

    identity = matches / aligned
    coverage = aligned / shorter
    return identity, coverage, aligned, shorter

def _verify_chunk(pairs: List[Tuple[int, int]]) -> List[Tuple[int, int, float, float, int, int]]:
    assert G_SEQS is not None
    out = []
    for i, j in pairs:
        identity, coverage, aligned_res, len_shorter = _nw_stats(G_SEQS[i], G_SEQS[j])
        if identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE:
            out.append((i, j, identity, coverage, aligned_res, len_shorter))
    return out

def _chunk(items: List[Tuple[int, int]], size: int):
    for k in range(0, len(items), size):
        yield items[k:k+size]

def main():
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f"Missing FASTA: {FASTA_PATH}")
    if not MMSEQS_TSV.exists():
        raise FileNotFoundError(f"Missing MMseqs TSV: {MMSEQS_TSV}")

    # FASTA numeric IDs by order
    records = list(SeqIO.parse(str(FASTA_PATH), "fasta"))
    ids = [r.id for r in records]
    seqs = [str(r.seq) for r in records]
    lens = [len(s) for s in seqs]
    id_to_idx: Dict[str, int] = {fid: i for i, fid in enumerate(ids)}

    # Read edges from mmseqs rep->member
    edges_set = set()
    skipped = 0
    with MMSEQS_TSV.open("r", encoding="utf-8") as f:
        for line in f:
            rep, mem = line.rstrip("\n").split("\t")
            if rep not in id_to_idx or mem not in id_to_idx:
                skipped += 1
                continue
            i = id_to_idx[rep]
            j = id_to_idx[mem]
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            edges_set.add((a, b))

    edges = sorted(edges_set)
    print("Total chains:", len(ids))
    print("MMseqs candidate edges:", len(edges))
    if skipped:
        print("Skipped mmseqs lines (IDs not found):", skipped)

    out_path = OUT_DIR / "nw_verified_pairs.tsv"
    tmp_path = OUT_DIR / "nw_verified_pairs.tsv.tmp"

    with tmp_path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("i\tj\tid_i\tid_j\tidentity\tcoverage_shorter\taligned_res\tlen_shorter\tlen_i\tlen_j\n")

    t0 = time.time()
    ctx = mp.get_context("spawn")

    passed_total = 0
    with ctx.Pool(
        processes=NW_WORKERS,
        initializer=_init_worker,
        initargs=(seqs, ids, lens, GAP_OPEN, GAP_EXTEND),
        maxtasksperchild=300,
    ) as pool:
        for verified in pool.imap_unordered(_verify_chunk, _chunk(edges, PAIRS_PER_CHUNK), chunksize=1):
            if not verified:
                continue
            with tmp_path.open("a", encoding="utf-8", newline="\n") as out:
                for (i, j, identity, coverage, aligned_res, len_shorter) in verified:
                    out.write(
                        f"{i}\t{j}\t{ids[i]}\t{ids[j]}\t"
                        f"{identity:.6f}\t{coverage:.6f}\t{aligned_res}\t{len_shorter}\t{lens[i]}\t{lens[j]}\n"
                    )
                    passed_total += 1

    if out_path.exists():
        out_path.unlink()
    tmp_path.rename(out_path)

    elapsed = time.time() - t0
    print("\n===== STEP 3.1–3.2 DONE =====")
    print("NW-verified edges:", passed_total)
    print("Wrote:", out_path)
    print("Elapsed:", f"{elapsed/3600:.2f} hours")

if __name__ == "__main__":
    main()
