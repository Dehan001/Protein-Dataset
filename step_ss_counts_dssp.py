from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import os
import threading
import time
from datetime import datetime

from Bio import SeqIO

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"
DSSP_DIR = BASE / "dssp_out"
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

ALL_CHAINS_FASTA = OUT_DIR / "all_chains.fasta"

MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")
LIBCIFPP_DATA_DIR = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\share\libcifpp")

WORKERS = 8
MKDSSP_TIMEOUT_SEC = 900

# DSSP codes
BETA_SS = {"E", "B"}          # beta strand + beta bridge
HELIX_SS = {"H", "G", "I"}    # alpha, 3-10, pi helix

# NEW outputs (no overwrite)
RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S")
SS_COUNTS_TSV = OUT_DIR / f"ss_counts_{RUN_ID}.tsv"
FAIL_LOG = OUT_DIR / f"dssp_failures_{RUN_ID}.log"
# ====================================================


def build_fasta_id_map(fasta_path: Path) -> Dict[Tuple[str, str], list[str]]:
    """
    Map (PDB_ID, CHAIN) -> [FASTA_IDs]
    Assumes FASTA IDs like: PDBID_..._CHAIN (chain = last token after underscore)
    """
    m: Dict[Tuple[str, str], list[str]] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        fid = rec.id.strip()
        parts = fid.split("_")
        if len(parts) < 2:
            continue
        pdb = parts[0].upper()
        chain = parts[-1]
        m.setdefault((pdb, chain), []).append(fid)
    return m


def fasta_ids_by_chain_for_pdb(fasta_map: Dict[Tuple[str, str], list[str]], pdb_id: str) -> Dict[str, list[str]]:
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
    Try DSSP text first; if mkdssp requests mmCIF output, rerun in mmCIF mode.
    """
    env = os.environ.copy()
    env["LIBCIFPP_DATA_DIR"] = str(LIBCIFPP_DATA_DIR)

    rc, so, se, _ = run_cmd([str(MKDSSP_EXE), str(cif_path), str(out_text_path)], env, MKDSSP_TIMEOUT_SEC)
    if rc == 0 and out_text_path.exists():
        return True, out_text_path, "dssp", ""

    err = (se or "") + (so or "")

    if "old DSSP format" in err or "mmCIF format instead" in err:
        rc2, so2, se2, _ = run_cmd(
            [str(MKDSSP_EXE), "--output-format=mmcif", str(cif_path), str(out_cif_path)],
            env,
            MKDSSP_TIMEOUT_SEC,
        )
        if rc2 == 0 and out_cif_path.exists():
            return True, out_cif_path, "mmcif", ""
        return False, out_cif_path, "mmcif", (se2 or "") + (so2 or "")

    return False, out_text_path, "dssp", err


def count_from_dssp_text(dssp_path: Path) -> Dict[str, Tuple[int, int, int, int]]:
    """
    Return per-chain:
      beta_res, helix_res, beta_strands, helices
    """
    lines = dssp_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    start = None
    for idx, line in enumerate(lines):
        if line.startswith("  #  RESIDUE"):
            start = idx + 1
            break
    if start is None:
        return {}

    beta_res: Dict[str, int] = {}
    helix_res: Dict[str, int] = {}
    beta_strands: Dict[str, int] = {}
    helices: Dict[str, int] = {}

    prev_beta: Dict[str, bool] = {}
    prev_helix: Dict[str, bool] = {}

    for line in lines[start:]:
        if len(line) < 17:
            continue
        chain = line[11].strip()
        ss = line[16].strip()
        if not chain:
            continue

        # Residue counts
        if ss in BETA_SS:
            beta_res[chain] = beta_res.get(chain, 0) + 1
        if ss in HELIX_SS:
            helix_res[chain] = helix_res.get(chain, 0) + 1

        # Segment counts
        is_beta = ss in BETA_SS
        was_beta = prev_beta.get(chain, False)
        if is_beta and not was_beta:
            beta_strands[chain] = beta_strands.get(chain, 0) + 1
        prev_beta[chain] = is_beta

        is_helix = ss in HELIX_SS
        was_helix = prev_helix.get(chain, False)
        if is_helix and not was_helix:
            helices[chain] = helices.get(chain, 0) + 1
        prev_helix[chain] = is_helix

    out: Dict[str, Tuple[int, int, int, int]] = {}
    chains = set(beta_res) | set(helix_res) | set(beta_strands) | set(helices)
    for ch in chains:
        out[ch] = (
            beta_res.get(ch, 0),
            helix_res.get(ch, 0),
            beta_strands.get(ch, 0),
            helices.get(ch, 0),
        )
    return out


def process_one(cif_path: Path) -> tuple[str, bool, str, Dict[str, Tuple[int,int,int,int]], str]:
    pdb_id = cif_path.stem.upper()
    out_text = DSSP_DIR / f"{pdb_id}.dssp"
    out_cif = DSSP_DIR / f"{pdb_id}.dssp.cif"
    ok, out_used, mode, err = run_mkdssp_adaptive(cif_path, out_text, out_cif)
    if not ok:
        return pdb_id, False, mode, {}, err
    if mode == "dssp":
        counts = count_from_dssp_text(out_used)
    else:
        # If you want, we can extend mmCIF parsing later.
        # For now: if mmCIF mode happens, we keep empty -> downstream will default zeros.
        counts = {}
    return pdb_id, True, mode, counts, ""


def main():
    if not ALL_CHAINS_FASTA.exists():
        raise FileNotFoundError(ALL_CHAINS_FASTA)
    if not MKDSSP_EXE.exists():
        raise FileNotFoundError(MKDSSP_EXE)
    if not LIBCIFPP_DATA_DIR.exists():
        raise FileNotFoundError(LIBCIFPP_DATA_DIR)

    cif_files = sorted(PDB_MODELS_DIR.glob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found in {PDB_MODELS_DIR}")

    fasta_map = build_fasta_id_map(ALL_CHAINS_FASTA)

    FAIL_LOG.write_text("", encoding="utf-8")
    SS_COUNTS_TSV.write_text(
        "fasta_id\tbeta_res\thelix_res\tbeta_strands\thelices\n",
        encoding="utf-8"
    )

    tsv_lock = threading.Lock()
    start_time = time.time()
    ok = fail = 0

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process_one, cif): cif for cif in cif_files}
        done = 0

        for fut in as_completed(futs):
            done += 1
            pdb_id, success, mode, counts_by_chain, err = fut.result()
            chain_to_fasta_ids = fasta_ids_by_chain_for_pdb(fasta_map, pdb_id)

            if not chain_to_fasta_ids:
                if not success:
                    fail += 1
                    with FAIL_LOG.open("a", encoding="utf-8") as f:
                        f.write(f"[FAIL] {pdb_id} mode={mode} (no FASTA mapping)\n{err}\n{'-'*80}\n")
                else:
                    ok += 1
                continue

            if not success:
                fail += 1
                with FAIL_LOG.open("a", encoding="utf-8") as f:
                    f.write(f"[FAIL] {pdb_id} mode={mode}\n{err}\n{'-'*80}\n")
                # default zeros for all mapped FASTA IDs
                counts_by_chain = {}

            else:
                ok += 1

            rows = []
            for ch, fids in chain_to_fasta_ids.items():
                beta_res, helix_res, beta_strands, helices = counts_by_chain.get(ch, (0, 0, 0, 0))
                for fid in fids:
                    rows.append((fid, beta_res, helix_res, beta_strands, helices))

            with tsv_lock:
                with SS_COUNTS_TSV.open("a", encoding="utf-8", newline="\n") as out:
                    for fid, br, hr, bs, hx in sorted(rows):
                        out.write(f"{fid}\t{br}\t{hr}\t{bs}\t{hx}\n")

            if done % 20 == 0 or done == len(cif_files):
                elapsed = time.time() - start_time
                print(f"Progress: {done}/{len(cif_files)} ok={ok} fail={fail} elapsed={elapsed:.1f}s")

    print("\nDONE")
    print("Wrote:", SS_COUNTS_TSV)
    print("Fail log:", FAIL_LOG)


if __name__ == "__main__":
    main()

# .\venv\Scripts\python.exe .\step_ss_counts_dssp.py 