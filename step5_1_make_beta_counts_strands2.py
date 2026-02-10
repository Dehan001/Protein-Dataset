from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import os
import threading
import time

from Bio import SeqIO  # pip install biopython

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"         # your *.cif coordinate files
DSSP_DIR = BASE / "dssp_out"                 # where we store mkdssp outputs
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

# FASTA with IDs like: PDBID_EMD-xxxxx_A
ALL_CHAINS_FASTA = OUT_DIR / "all_chains.fasta"

# mkdssp that you confirmed works when called by full path
MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")

# Important for the "mmcif_pdbx.dic" issue
LIBCIFPP_DATA_DIR = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\share\libcifpp")

WORKERS = 8
MKDSSP_TIMEOUT_SEC = 900  # 15 min per file

BETA_COUNTS_TSV = OUT_DIR / "beta_counts.tsv"
FAIL_LOG = OUT_DIR / "dssp_failures.log"
ZEROCOUNT_LOG = OUT_DIR / "dssp_zero_counts.log"

# DSSP beta codes
BETA_SS = {"E", "B"}
# ====================================================


def build_fasta_id_map(fasta_path: Path) -> Dict[Tuple[str, str], list[str]]:
    """
    Map (PDB_ID, CHAIN) -> [FASTA_IDs]
    FASTA IDs expected: PDBID_EMD-xxxxx_CHAIN
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing {fasta_path}")

    m: Dict[Tuple[str, str], list[str]] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        fid = rec.id.strip()
        parts = fid.split("_")
        if len(parts) < 2:
            continue
        pdb = parts[0].upper()
        chain = parts[-1]  # last token is chain ID
        key = (pdb, chain)
        m.setdefault(key, []).append(fid)
    return m


def fasta_ids_by_chain_for_pdb(
    fasta_map: Dict[Tuple[str, str], list[str]],
    pdb_id: str
) -> Dict[str, list[str]]:
    """Return chain -> [fasta_ids] for a given PDB."""
    pdb_id = pdb_id.upper()
    out: Dict[str, list[str]] = {}
    for (pdb, chain), ids in fasta_map.items():
        if pdb == pdb_id:
            out.setdefault(chain, []).extend(ids)
    return out


def run_cmd(cmd: list[str], env: dict, timeout_sec: int) -> tuple[int, str, str, bool]:
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    try:
        so, se = p.communicate(timeout=timeout_sec)
        return p.returncode, so, se, False
    except subprocess.TimeoutExpired:
        try:
            p.kill()
        except Exception:
            pass
        so, se = p.communicate()
        return 124, so or "", (se or "") + f"\n[TIMEOUT after {timeout_sec}s]", True


def run_mkdssp_adaptive(cif_path: Path, out_text_path: Path, out_cif_path: Path) -> tuple[bool, Path, str, str]:
    """
    Try classic DSSP text first; if mkdssp says "use mmCIF format", rerun to produce mmCIF output.
    Returns: (ok, output_path_used, mode, err_msg) where mode in {"dssp","mmcif"}
    """
    env = os.environ.copy()
    env["LIBCIFPP_DATA_DIR"] = str(LIBCIFPP_DATA_DIR)

    # 1) classic DSSP text: mkdssp input output
    rc, so, se, _ = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_text_path)], env, MKDSSP_TIMEOUT_SEC)
    if rc == 0 and out_text_path.exists():
        return True, out_text_path, "dssp", ""

    err = (se or "") + (so or "")

    # 2) switch to mmCIF output if needed
    if "old DSSP format" in err or "mmCIF format instead" in err:
        rc2, so2, se2, _ = run_cmd(
            [str(MKDSSP_EXE), "--output-format=mmcif", str(cif_path), str(out_cif_path)],
            env,
            MKDSSP_TIMEOUT_SEC,
        )
        if rc2 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif", ""

        # fallback: infer from extension
        rc3, so3, se3, _ = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_cif_path)], env, MKDSSP_TIMEOUT_SEC)
        if rc3 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif", ""

        return False, out_cif_path, "mmcif", (se2 or "") + (so2 or "") + (se3 or "") + (so3 or "")

    return False, out_text_path, "dssp", err


def count_beta_strands_from_dssp_text(dssp_path: Path) -> Dict[str, int]:
    """
    Count β-STRANDS (segments), not residues.
    For each chain: increment when SS enters {E,B} from non-{E,B}.
    """
    lines = dssp_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    start = None
    for idx, line in enumerate(lines):
        if line.startswith("  #  RESIDUE"):
            start = idx + 1
            break
    if start is None:
        return {}

    beta_strands: Dict[str, int] = {}
    prev_is_beta: Dict[str, bool] = {}

    for line in lines[start:]:
        if len(line) < 17:
            continue
        chain = line[11].strip()
        ss = line[16].strip()
        if not chain:
            continue

        is_beta = (ss in BETA_SS)
        was_beta = prev_is_beta.get(chain, False)
        if is_beta and not was_beta:
            beta_strands[chain] = beta_strands.get(chain, 0) + 1
        prev_is_beta[chain] = is_beta

    return beta_strands


# ---------------- mmCIF text loop parser (no gemmi needed) ----------------

def _is_mmcif_tag_line(line: str) -> bool:
    return line.startswith("_")


def _tokenize_mmcif_row(line: str) -> List[str]:
    """
    Simple tokenizer for mkdssp mmCIF: mkdssp outputs are typically whitespace-separated.
    We avoid full mmCIF quoting complexity; if a row has quotes, split() still often works for mkdssp output.
    """
    return line.strip().split()


def parse_mkdssp_mmcif_beta_strands(mmcif_path: Path) -> Dict[str, int]:
    """
    Parse mkdssp mmCIF output and count β-strand segments per chain (transitions into E/B).

    Strategy:
    - Scan loop_ blocks.
    - Collect tags.
    - If tags include something like dssp/sec_struct/secondary_structure, try to locate:
        * chain column: auth_asym_id / label_asym_id / asym_id / chain_id
        * ss column: dssp / sec_struct / secondary_structure / ss
    - Then count transitions into E/B per chain.
    """
    try:
        text = mmcif_path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return {}

    lines = [ln.rstrip() for ln in text.splitlines()]

    i = 0
    n = len(lines)
    best: Dict[str, int] = {}

    while i < n:
        line = lines[i].strip()
        if line.lower() != "loop_":
            i += 1
            continue

        # Collect tags
        i += 1
        tags: List[str] = []
        while i < n:
            ln = lines[i].strip()
            if not ln:
                i += 1
                continue
            if _is_mmcif_tag_line(ln):
                tags.append(ln)
                i += 1
                continue
            break  # first data row or something else

        if not tags:
            continue

        tags_l = [t.lower() for t in tags]

        # must "smell" like DSSP/sec-struct
        if not any(("dssp" in t or "sec_struct" in t or "secondary_structure" in t) for t in tags_l):
            # skip data rows until loop ends
            while i < n:
                ln = lines[i].strip()
                if not ln or ln.lower() == "loop_" or _is_mmcif_tag_line(ln) or ln.startswith("#"):
                    break
                i += 1
            continue

        # find chain column
        chain_idx: Optional[int] = None
        for idx, t in enumerate(tags_l):
            if any(k in t for k in ("auth_asym_id", "label_asym_id", "asym_id", "chain_id")):
                chain_idx = idx
                break

        # find ss column
        ss_idx: Optional[int] = None
        for idx, t in enumerate(tags_l):
            if any(k in t for k in ("secondary_structure", "sec_struct", ".ss", "dssp")):
                ss_idx = idx
                break

        if chain_idx is None or ss_idx is None:
            # skip data rows until loop ends
            while i < n:
                ln = lines[i].strip()
                if not ln or ln.lower() == "loop_" or _is_mmcif_tag_line(ln) or ln.startswith("#"):
                    break
                i += 1
            continue

        beta_strands: Dict[str, int] = {}
        prev_is_beta: Dict[str, bool] = {}

        # Read data rows
        while i < n:
            ln = lines[i].strip()
            if not ln:
                i += 1
                continue
            if ln.lower() == "loop_" or _is_mmcif_tag_line(ln) or ln.startswith("#"):
                break

            row = _tokenize_mmcif_row(ln)
            # If row is shorter than expected, skip
            if len(row) <= max(chain_idx, ss_idx):
                i += 1
                continue

            chain = str(row[chain_idx]).strip().strip('"').strip("'")
            ss = str(row[ss_idx]).strip().strip('"').strip("'")

            if chain and ss and ss not in (".", "?"):
                code = ss[0]
                is_beta = (code in BETA_SS)
                was_beta = prev_is_beta.get(chain, False)
                if is_beta and not was_beta:
                    beta_strands[chain] = beta_strands.get(chain, 0) + 1
                prev_is_beta[chain] = is_beta

            i += 1

        if beta_strands:
            return beta_strands

        best = beta_strands  # keep last attempt

    return best


def process_one(cif_path: Path) -> tuple[str, bool, str, Dict[str, int], str]:
    """
    Returns: (pdb_id, success, mode, counts_by_chain, err)
    """
    pdb_id = cif_path.stem.upper()
    out_text = DSSP_DIR / f"{pdb_id}.dssp"
    out_cif = DSSP_DIR / f"{pdb_id}.dssp.cif"

    try:
        ok, out_used, mode, err = run_mkdssp_adaptive(cif_path, out_text, out_cif)
        if not ok:
            return pdb_id, False, mode, {}, err

        if mode == "dssp":
            counts = count_beta_strands_from_dssp_text(out_used)
        else:
            counts = parse_mkdssp_mmcif_beta_strands(out_used)

        # If mmCIF parsing fails silently (empty), still "success"—main will assign zeros.
        return pdb_id, True, mode, counts, ""
    except Exception as e:
        # Never crash the whole pool
        return pdb_id, False, "exception", {}, f"[EXCEPTION] {type(e).__name__}: {e}"


def main():
    cif_files = sorted(PDB_MODELS_DIR.glob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found in {PDB_MODELS_DIR}")
    if not MKDSSP_EXE.exists():
        raise FileNotFoundError(f"mkdssp.exe not found: {MKDSSP_EXE}")
    if not LIBCIFPP_DATA_DIR.exists():
        raise FileNotFoundError(f"LIBCIFPP_DATA_DIR not found: {LIBCIFPP_DATA_DIR}")

    fasta_map = build_fasta_id_map(ALL_CHAINS_FASTA)

    print("===== STEP 5.1 (DSSP -> beta STRANDS per FASTA chain) =====")
    print("Total CIFs:", len(cif_files))
    print("Workers:", WORKERS)

    FAIL_LOG.write_text("", encoding="utf-8")
    ZEROCOUNT_LOG.write_text("", encoding="utf-8")
    BETA_COUNTS_TSV.write_text("fasta_id\tbeta_strands\n", encoding="utf-8")
    tsv_lock = threading.Lock()

    ok = 0
    fail = 0
    zero = 0
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process_one, cif): cif for cif in cif_files}

        done = 0
        for fut in as_completed(futs):
            done += 1
            pdb_id, success, mode, counts_by_chain, err = fut.result()

            # All FASTA IDs that belong to this PDB (from all_chains.fasta)
            chain_to_fasta_ids = fasta_ids_by_chain_for_pdb(fasta_map, pdb_id)

            # If this PDB doesn't exist in the FASTA map, there's nothing to write
            if not chain_to_fasta_ids:
                if not success:
                    fail += 1
                    with FAIL_LOG.open("a", encoding="utf-8") as f:
                        f.write(f"[FAIL] {pdb_id} mode={mode} (no FASTA IDs found for this PDB)\n{err}\n{'-'*80}\n")
                else:
                    ok += 1
                continue

            # Default: assign 0 beta strands to every chain in FASTA
            final_counts_by_chain: Dict[str, int] = {ch: 0 for ch in chain_to_fasta_ids.keys()}

            if not success:
                fail += 1
                # keep all zeros; just log the failure
                with FAIL_LOG.open("a", encoding="utf-8") as f:
                    f.write(f"[FAIL] {pdb_id} mode={mode}\n{err}\n{'-'*80}\n")
            else:
                ok += 1
                # overwrite zeros with actual counts where available
                for ch, b in counts_by_chain.items():
                    if ch in final_counts_by_chain:
                        final_counts_by_chain[ch] = b

            # If everything is still zero, log it (covers: true zero beta, mmcif-parse-empty, and failures)
            if all(v == 0 for v in final_counts_by_chain.values()):
                zero += 1
                with ZEROCOUNT_LOG.open("a", encoding="utf-8") as f:
                    f.write(f"[ZERO] {pdb_id} mode={mode} (all FASTA chains assigned 0)\n")

            # Write a row for EVERY FASTA ID for this PDB (chain-level), always.
            rows = []
            for ch, fasta_ids in chain_to_fasta_ids.items():
                b = final_counts_by_chain.get(ch, 0)
                for fid in fasta_ids:
                    rows.append((fid, b))

            with tsv_lock:
                with BETA_COUNTS_TSV.open("a", encoding="utf-8", newline="\n") as out:
                    for fid, b in sorted(rows):
                        out.write(f"{fid}\t{b}\n")

            if done % 10 == 0 or done == len(cif_files):
                elapsed = time.time() - start_time
                print(f"Progress: {done}/{len(cif_files)} (ok={ok}, fail={fail}, zero={zero}) elapsed={elapsed:.1f}s")
    
    # BETA_COUNTS_TSV.write_text("fasta_id\tbeta_strands\n", encoding="utf-8")

    print("\n===== STEP 5.1 DONE =====")
    print("OK:", ok, "FAIL:", fail, "ZERO(all chains 0):", zero)
    print("Wrote:", BETA_COUNTS_TSV)
    print("Fail log:", FAIL_LOG)
    print("Zero log:", ZEROCOUNT_LOG)


if __name__ == "__main__":
    main()

# in your project folder
# .\venv\Scripts\python.exe -m pip install biopython
# (gemmi no longer required in this version)
# .\venv\Scripts\python.exe .\step5_1_make_beta_counts_strands2.py
