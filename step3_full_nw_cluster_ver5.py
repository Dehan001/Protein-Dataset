from __future__ import annotations

from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
import multiprocessing as mp
import time
import os

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
FASTA_PATH = BASE / "outputs" / "all_chains.fasta"
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)

# Similarity criteria (Step 2)
MIN_IDENTITY = 0.90
MIN_COVERAGE = 0.80

# NW config (Step 1)
GAP_OPEN = -10.0
GAP_EXTEND = -0.5

# Workers (your supervisor: 8–10 threads/processes)
NW_WORKERS = 8  # set 10 if you want

# Chunk size: number of (i,j) pairs per worker task
# Bigger = fewer overhead calls, but more RAM usage per task
PAIRS_PER_CHUNK = 2500

# Optional fast filter to reduce compute (SAFE with coverage>=0.8):
# If max_len/min_len > 1/0.8 = 1.25, then even perfect alignment can't reach 0.8 coverage of shorter
# This is mathematically safe and does NOT change correctness for your thresholds.
USE_LENGTH_RATIO_FILTER = True
LENGTH_RATIO_MAX = 1.25  # derived from 0.8 coverage

# Progress printing
PRINT_EVERY_CHUNKS = 80
# ====================================================

# ------------ Worker globals ------------
G_SEQS: List[str] | None = None
G_LENS: List[int] | None = None
G_BLOSUM62 = None
G_GAP_OPEN: float | None = None
G_GAP_EXTEND: float | None = None

def _init_worker(seqs: List[str], lens: List[int], gap_open: float, gap_extend: float):
    global G_SEQS, G_LENS, G_BLOSUM62, G_GAP_OPEN, G_GAP_EXTEND
    G_SEQS = seqs
    G_LENS = lens
    G_GAP_OPEN = gap_open
    G_GAP_EXTEND = gap_extend
    # Load once per worker
    G_BLOSUM62 = substitution_matrices.load("BLOSUM62")

def _nw_identity_coverage(seq1: str, seq2: str) -> Tuple[float, float, int, int]:
    """
    Needleman-Wunsch global alignment (pairwise2.globalds) with BLOSUM62.

    identity = matches / aligned_residues
    coverage_shorter = aligned_residues / len(shorter_chain)
    aligned_residues = positions where both residues are not gaps
    """
    alns = pairwise2.align.globalds(
        seq1, seq2, G_BLOSUM62, G_GAP_OPEN, G_GAP_EXTEND, one_alignment_only=True
    )
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

def _process_chunk(pairs: List[Tuple[int, int]]) -> List[Tuple[int, int, float, float, int, int]]:
    """
    Worker: compute NW stats for a chunk of pairs, return only passing edges.
    Returns tuples: (i, j, identity, coverage, aligned_res, len_shorter)
    """
    assert G_SEQS is not None
    out = []
    for i, j in pairs:
        identity, coverage, aligned_res, len_shorter = _nw_identity_coverage(G_SEQS[i], G_SEQS[j])
        if identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE:
            out.append((i, j, identity, coverage, aligned_res, len_shorter))
    return out

def _pair_chunk_generator(n: int, lens: List[int]):
    """
    Generate chunks of pairs (i,j) for all i<j.
    Applies SAFE length ratio filter if enabled.
    """
    chunk: List[Tuple[int, int]] = []
    tested_in_chunk = 0

    for i in range(n):
        li = lens[i]
        for j in range(i + 1, n):
            lj = lens[j]

            if USE_LENGTH_RATIO_FILTER:
                mn = li if li < lj else lj
                mx = lj if lj > li else li
                if mn == 0:
                    continue
                if (mx / mn) > LENGTH_RATIO_MAX:
                    continue

            chunk.append((i, j))
            tested_in_chunk += 1
            if len(chunk) >= PAIRS_PER_CHUNK:
                yield chunk
                chunk = []
                tested_in_chunk = 0

    if chunk:
        yield chunk

def main():
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f"Missing FASTA: {FASTA_PATH}")

    records = list(SeqIO.parse(str(FASTA_PATH), "fasta"))
    ids = [r.id for r in records]
    seqs = [str(r.seq) for r in records]
    lens = [len(s) for s in seqs]
    n = len(seqs)

    print("===== FULL ALL-PAIRS NW (BLOSUM62) =====")
    print("FASTA:", FASTA_PATH)
    print("Total sequences:", n)
    print(f"Workers: {NW_WORKERS}")
    print(f"Thresholds: identity>={MIN_IDENTITY}, coverage>={MIN_COVERAGE}")
    if USE_LENGTH_RATIO_FILTER:
        print(f"SAFE length-ratio filter ON: max_len/min_len <= {LENGTH_RATIO_MAX} (derived from coverage>=0.8)")
    else:
        print("Length-ratio filter OFF (full brute force)")

    # Output: ONLY passing edges
    out_path = OUT_DIR / "similar_pairs.tsv"
    tmp_path = OUT_DIR / "similar_pairs.tsv.tmp"

    with tmp_path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("i\tj\tid_i\tid_j\tidentity\tcoverage_shorter\taligned_res\tlen_shorter\tlen_i\tlen_j\n")

    ctx = mp.get_context("spawn")
    t0 = time.time()

    chunks_done = 0
    passing_edges = 0

    # Stream computation
    with ctx.Pool(
        processes=NW_WORKERS,
        initializer=_init_worker,
        initargs=(seqs, lens, GAP_OPEN, GAP_EXTEND),
        maxtasksperchild=300,
    ) as pool:

        for passed_list in pool.imap_unordered(_process_chunk, _pair_chunk_generator(n, lens), chunksize=1):
            chunks_done += 1

            if passed_list:
                with tmp_path.open("a", encoding="utf-8", newline="\n") as out:
                    for (i, j, identity, coverage, aligned_res, len_shorter) in passed_list:
                        out.write(
                            f"{i}\t{j}\t{ids[i]}\t{ids[j]}\t"
                            f"{identity:.6f}\t{coverage:.6f}\t{aligned_res}\t{len_shorter}\t{lens[i]}\t{lens[j]}\n"
                        )
                        passing_edges += 1

            if chunks_done % PRINT_EVERY_CHUNKS == 0:
                elapsed = time.time() - t0
                print(f"Chunks done: {chunks_done} | passing_edges={passing_edges} | elapsed={elapsed/3600:.2f}h")

    # Finalize output atomically
    if out_path.exists():
        out_path.unlink()
    tmp_path.rename(out_path)

    elapsed = time.time() - t0
    print("\n===== STEP 1–2 DONE (FULL) =====")
    print("Passing edges written:", passing_edges)
    print("Wrote:", out_path)
    print("Elapsed:", f"{elapsed/3600:.2f} hours")


if __name__ == "__main__":
    main()
