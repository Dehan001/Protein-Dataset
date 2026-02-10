from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, List, Optional, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import os
import threading
import time
from datetime import datetime

from Bio import SeqIO  # pip install biopython
import shlex

# ===================== SETTINGS =====================
BASE = Path(__file__).resolve().parent
PDB_MODELS_DIR = BASE / "pdb_models"     # structure mmCIF files
DSSP_DIR = BASE / "dssp_out"             # mkdssp outputs (text or mmcif)
OUT_DIR = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)
DSSP_DIR.mkdir(exist_ok=True)

ALL_CHAINS_FASTA = OUT_DIR / "all_chains.fasta"

# mkdssp
MKDSSP_EXE = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\bin\mkdssp.exe")
LIBCIFPP_DATA_DIR = Path(r"C:\Users\fd02629\AppData\Local\anaconda3\share\libcifpp")

WORKERS = 8
MKDSSP_TIMEOUT_SEC = 900   # you can lower to 300 if some are “stuck”

# DSSP codes
BETA_SS = {"E", "B"}          # beta strand + beta bridge
HELIX_SS = {"H", "G", "I"}    # alpha, 3-10, pi helix

# Output files (never overwrite: timestamped)
RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S")
SS_COUNTS_TSV = OUT_DIR / f"ss_counts_{RUN_ID}.tsv"
FAIL_LOG = OUT_DIR / f"dssp_failures_{RUN_ID}.log"
# ====================================================


# ----------------- FASTA mapping -----------------

def pdb_from_fasta_id(fid: str) -> str:
    return fid.split("_", 1)[0].upper()

def chain_from_fasta_id(fid: str) -> str:
    # your IDs look like PDBID_EMD-xxxxx_CHAIN
    return fid.rsplit("_", 1)[-1].strip()

def build_fasta_id_map(fasta_path: Path) -> Dict[str, Dict[str, List[str]]]:
    """
    Returns:
      pdb_id -> { chain_id_in_fasta -> [fasta_ids...] }
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing FASTA: {fasta_path}")

    out: Dict[str, Dict[str, List[str]]] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        fid = rec.id.strip()
        pdb = pdb_from_fasta_id(fid)
        ch = chain_from_fasta_id(fid)
        out.setdefault(pdb, {}).setdefault(ch, []).append(fid)
    return out


# ----------------- locate structure CIF per PDB -----------------

def index_cif_files(models_dir: Path) -> Dict[str, Path]:
    """
    Build mapping PDBID -> structure cif path by searching filename containing the PDBID.
    Works even if filenames are like: "8GHA_model_1.cif", "pdb_8gha_xxx.cif", etc.
    If multiple match, we pick the shortest filename (usually the clean one).
    """
    cifs = list(models_dir.glob("*.cif"))
    if not cifs:
        raise FileNotFoundError(f"No *.cif found in {models_dir}")

    # Precompute uppercase names
    name_to_path = [(p.name.upper(), p) for p in cifs]

    def find_for_pdb(pdb: str) -> Optional[Path]:
        pdb_u = pdb.upper()
        matches = [path for name_u, path in name_to_path if pdb_u in name_u]
        if not matches:
            return None
        matches.sort(key=lambda x: (len(x.name), x.name))
        return matches[0]

    # We return a function-like mapping by filling later per pdb
    # (We keep the list and do find_for_pdb on demand)
    # But returning a dict is convenient too:
    return { }  # placeholder; we'll use find_for_pdb closure below


# ----------------- mkdssp runner -----------------

def run_cmd(cmd: list[str], env: dict, timeout_sec: int) -> Tuple[int, str, str, bool]:
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


def run_mkdssp_adaptive(struct_cif: Path, out_text: Path, out_mmcif: Path) -> Tuple[bool, Path, str, str]:
    """
    Try classic DSSP text; if mkdssp says use mmCIF output, rerun with --output-format=mmcif.
    Returns: ok, output_path_used, mode("dssp"/"mmcif"/"timeout"), err
    """
    env = os.environ.copy()
    env["LIBCIFPP_DATA_DIR"] = str(LIBCIFPP_DATA_DIR)

    rc, so, se, timed = run_cmd([str(MKDSSP_EXE), str(struct_cif), str(out_text)], env, MKDSSP_TIMEOUT_SEC)
    if rc == 0 and out_text.exists():
        return True, out_text, "dssp", ""
    if timed:
        return False, out_text, "timeout", (se or "") + (so or "")

    err = (se or "") + (so or "")
    if ("old DSSP format" in err) or ("mmCIF format instead" in err) or ("--output-format=mmcif" in err):
        rc2, so2, se2, timed2 = run_cmd(
            [str(MKDSSP_EXE), "--output-format=mmcif", str(struct_cif), str(out_mmcif)],
            env,
            MKDSSP_TIMEOUT_SEC
        )
        if rc2 == 0 and out_mmcif.exists():
            return True, out_mmcif, "mmcif", ""
        if timed2:
            return False, out_mmcif, "timeout", (se2 or "") + (so2 or "")
        return False, out_mmcif, "mmcif", (se2 or "") + (so2 or "")

    return False, out_text, "dssp", err


# ----------------- parsing DSSP text output -----------------

def parse_dssp_text_counts(dssp_path: Path) -> Dict[str, Tuple[int,int,int,int]]:
    """
    Returns per DSSP-chain-id:
      beta_res, helix_res, beta_strands, helices
    """
    lines = dssp_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("  #  RESIDUE"):
            start = i + 1
            break
    if start is None:
        return {}

    beta_res: Dict[str, int] = {}
    helix_res: Dict[str, int] = {}
    beta_strands: Dict[str, int] = {}
    helices: Dict[str, int] = {}
    prev_beta: Dict[str, bool] = {}
    prev_helix: Dict[str, bool] = {}

    for ln in lines[start:]:
        if len(ln) < 17:
            continue
        ch = ln[11].strip()
        ss = ln[16].strip()
        if not ch:
            continue

        # residue counts
        if ss in BETA_SS:
            beta_res[ch] = beta_res.get(ch, 0) + 1
        if ss in HELIX_SS:
            helix_res[ch] = helix_res.get(ch, 0) + 1

        # segment counts (transitions)
        is_beta = ss in BETA_SS
        was_beta = prev_beta.get(ch, False)
        if is_beta and not was_beta:
            beta_strands[ch] = beta_strands.get(ch, 0) + 1
        prev_beta[ch] = is_beta

        is_helix = ss in HELIX_SS
        was_helix = prev_helix.get(ch, False)
        if is_helix and not was_helix:
            helices[ch] = helices.get(ch, 0) + 1
        prev_helix[ch] = is_helix

    out: Dict[str, Tuple[int,int,int,int]] = {}
    chains = set(beta_res) | set(helix_res) | set(beta_strands) | set(helices)
    for ch in chains:
        out[ch] = (
            beta_res.get(ch, 0),
            helix_res.get(ch, 0),
            beta_strands.get(ch, 0),
            helices.get(ch, 0),
        )
    return out


# ----------------- minimal CIF loop parser helpers -----------------

def _tokenize_cif_row(line: str) -> List[str]:
    """
    Handles quotes reasonably well for typical mmCIF loop rows.
    """
    # shlex.split handles quoted strings
    try:
        return shlex.split(line, posix=True)
    except ValueError:
        # fallback
        return line.strip().split()

def _read_loop_table(lines: List[str], start_i: int) -> Tuple[List[str], List[List[str]], int]:
    """
    Given lines and index at 'loop_', read tags + rows.
    Returns tags, rows, next_index.
    """
    i = start_i + 1
    tags: List[str] = []
    n = len(lines)

    while i < n:
        ln = lines[i].strip()
        if not ln:
            i += 1
            continue
        if ln.startswith("_"):
            tags.append(ln)
            i += 1
            continue
        break

    rows: List[List[str]] = []
    while i < n:
        ln = lines[i].strip()
        if not ln:
            i += 1
            continue
        if ln.startswith("loop_") or ln.startswith("_") or ln.startswith("#"):
            break
        rows.append(_tokenize_cif_row(ln))
        i += 1

    return tags, rows, i


# ----------------- parse structure mmCIF for auth<->label chain mapping -----------------

def parse_auth_label_chain_map(struct_cif: Path) -> Tuple[Dict[str,str], Dict[str,str]]:
    """
    Parse structure mmCIF (_atom_site loop) to build:
      auth_asym_id -> label_asym_id
      label_asym_id -> auth_asym_id
    Uses the first observed mapping per ID.
    """
    try:
        text = struct_cif.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return {}, {}

    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    i = 0
    n = len(lines)

    a2l: Dict[str, str] = {}
    l2a: Dict[str, str] = {}

    while i < n:
        if lines[i].strip().lower() != "loop_":
            i += 1
            continue

        tags, rows, j = _read_loop_table(lines, i)
        tags_l = [t.lower() for t in tags]

        # Must be atom_site and contain both columns
        try:
            ia = tags_l.index("_atom_site.auth_asym_id")
            il = tags_l.index("_atom_site.label_asym_id")
        except ValueError:
            i = j
            continue

        for r in rows:
            if len(r) <= max(ia, il):
                continue
            auth = str(r[ia]).strip().strip('"').strip("'")
            lab = str(r[il]).strip().strip('"').strip("'")
            if auth and lab:
                if auth not in a2l:
                    a2l[auth] = lab
                if lab not in l2a:
                    l2a[lab] = auth

        # we can stop after we found a decent mapping
        if a2l:
            break

        i = j

    return a2l, l2a


# ----------------- parse mkdssp mmCIF output -----------------

def parse_dssp_mmcif_counts(dssp_mmcif: Path) -> Dict[str, Tuple[int,int,int,int]]:
    """
    Parse mkdssp mmCIF output. We look for a loop containing:
      - chain id column (auth_asym_id / label_asym_id / asym_id / chain_id)
      - ss column (dssp / secondary_structure / sec_struct / ss)
    Then compute residue and segment counts per chain.
    """
    try:
        text = dssp_mmcif.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return {}

    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    i = 0
    n = len(lines)

    while i < n:
        if lines[i].strip().lower() != "loop_":
            i += 1
            continue

        tags, rows, j = _read_loop_table(lines, i)
        tags_l = [t.lower() for t in tags]

        # only consider loops that “smell” like DSSP/sec-struct
        if not any(("dssp" in t or "sec_struct" in t or "secondary_structure" in t) for t in tags_l):
            i = j
            continue

        # chain col
        chain_idx: Optional[int] = None
        for k in ("auth_asym_id", "label_asym_id", "asym_id", "chain_id"):
            for idx, t in enumerate(tags_l):
                if k in t:
                    chain_idx = idx
                    break
            if chain_idx is not None:
                break

        # ss col
        ss_idx: Optional[int] = None
        for idx, t in enumerate(tags_l):
            if ("secondary_structure" in t) or ("sec_struct" in t) or t.endswith(".ss") or ("dssp" in t):
                ss_idx = idx
                break

        if chain_idx is None or ss_idx is None:
            i = j
            continue

        # Now compute counts
        beta_res: Dict[str, int] = {}
        helix_res: Dict[str, int] = {}
        beta_strands: Dict[str, int] = {}
        helices: Dict[str, int] = {}
        prev_beta: Dict[str, bool] = {}
        prev_helix: Dict[str, bool] = {}

        for r in rows:
            if len(r) <= max(chain_idx, ss_idx):
                continue
            ch = str(r[chain_idx]).strip().strip('"').strip("'")
            ss = str(r[ss_idx]).strip().strip('"').strip("'")
            if not ch or not ss or ss in (".", "?"):
                continue

            code = ss[0]

            # residue counts
            if code in BETA_SS:
                beta_res[ch] = beta_res.get(ch, 0) + 1
            if code in HELIX_SS:
                helix_res[ch] = helix_res.get(ch, 0) + 1

            # segment counts (transition)
            is_beta = code in BETA_SS
            was_beta = prev_beta.get(ch, False)
            if is_beta and not was_beta:
                beta_strands[ch] = beta_strands.get(ch, 0) + 1
            prev_beta[ch] = is_beta

            is_helix = code in HELIX_SS
            was_helix = prev_helix.get(ch, False)
            if is_helix and not was_helix:
                helices[ch] = helices.get(ch, 0) + 1
            prev_helix[ch] = is_helix

        out: Dict[str, Tuple[int,int,int,int]] = {}
        chains = set(beta_res) | set(helix_res) | set(beta_strands) | set(helices)
        if chains:
            for ch in chains:
                out[ch] = (
                    beta_res.get(ch, 0),
                    helix_res.get(ch, 0),
                    beta_strands.get(ch, 0),
                    helices.get(ch, 0),
                )
            return out

        i = j

    return {}


# ----------------- chain resolution (FASTA chain -> DSSP chain) -----------------

def resolve_chain_counts_for_fasta_chain(
    fasta_chain: str,
    counts_by_dssp_chain: Dict[str, Tuple[int,int,int,int]],
    auth2label: Dict[str, str],
    label2auth: Dict[str, str]
) -> Tuple[int,int,int,int]:
    """
    Try to match the FASTA chain to the DSSP chain keys.
    Order:
      1) direct (FASTA chain equals DSSP chain)
      2) map auth->label then lookup
      3) map label->auth then lookup
      4) 0s
    """
    if fasta_chain in counts_by_dssp_chain:
        return counts_by_dssp_chain[fasta_chain]

    lab = auth2label.get(fasta_chain)
    if lab and lab in counts_by_dssp_chain:
        return counts_by_dssp_chain[lab]

    auth = label2auth.get(fasta_chain)
    if auth and auth in counts_by_dssp_chain:
        return counts_by_dssp_chain[auth]

    return (0, 0, 0, 0)


# ----------------- worker per PDB -----------------

def process_one_pdb(pdb_id: str, struct_cif: Path, fasta_chains: Dict[str, List[str]]) -> Tuple[str, bool, str, List[Tuple[str,int,int,int,int]], str]:
    """
    Returns:
      pdb_id, success, mode, rows[(fasta_id, beta_res, helix_res, beta_strands, helices)], err
    """
    pdb_u = pdb_id.upper()
    out_text = DSSP_DIR / f"{pdb_u}.dssp"
    out_mmcif = DSSP_DIR / f"{pdb_u}.dssp.cif"

    ok, out_used, mode, err = run_mkdssp_adaptive(struct_cif, out_text, out_mmcif)
    if not ok:
        # still write zeros for this PDB's FASTA IDs
        rows: List[Tuple[str,int,int,int,int]] = []
        for ch, fids in fasta_chains.items():
            for fid in fids:
                rows.append((fid, 0, 0, 0, 0))
        return pdb_u, False, mode, rows, err

    # parse output
    if mode == "dssp":
        counts_by_chain = parse_dssp_text_counts(out_used)
    else:
        counts_by_chain = parse_dssp_mmcif_counts(out_used)

    # chain mapping from structure cif
    auth2label, label2auth = parse_auth_label_chain_map(struct_cif)

    # build output rows for every FASTA id
    rows: List[Tuple[str,int,int,int,int]] = []
    for fasta_chain, fids in fasta_chains.items():
        br, hr, bs, hx = resolve_chain_counts_for_fasta_chain(
            fasta_chain, counts_by_chain, auth2label, label2auth
        )
        for fid in fids:
            rows.append((fid, br, hr, bs, hx))

    return pdb_u, True, mode, rows, ""


# ----------------- main -----------------

def main():
    if not ALL_CHAINS_FASTA.exists():
        raise FileNotFoundError(ALL_CHAINS_FASTA)
    if not MKDSSP_EXE.exists():
        raise FileNotFoundError(MKDSSP_EXE)
    if not LIBCIFPP_DATA_DIR.exists():
        raise FileNotFoundError(LIBCIFPP_DATA_DIR)
    if not PDB_MODELS_DIR.exists():
        raise FileNotFoundError(PDB_MODELS_DIR)

    fasta_map = build_fasta_id_map(ALL_CHAINS_FASTA)  # pdb -> chain -> [fasta_ids]
    pdb_ids = sorted(fasta_map.keys())

    # index CIF filenames once (fast)
    all_cifs = list(PDB_MODELS_DIR.glob("*.cif"))
    name_to_path = [(p.name.upper(), p) for p in all_cifs]

    def find_struct_cif(pdb: str) -> Optional[Path]:
        pdb_u = pdb.upper()
        matches = [path for name_u, path in name_to_path if pdb_u in name_u]
        if not matches:
            return None
        matches.sort(key=lambda x: (len(x.name), x.name))
        return matches[0]

    # write headers (new file only; no overwrite)
    SS_COUNTS_TSV.write_text("fasta_id\tbeta_res\thelix_res\tbeta_strands\thelices\n", encoding="utf-8")
    FAIL_LOG.write_text("", encoding="utf-8")
    tsv_lock = threading.Lock()

    # Build jobs
    jobs: List[Tuple[str, Path, Dict[str, List[str]]]] = []
    missing_struct = 0
    for pdb in pdb_ids:
        cif = find_struct_cif(pdb)
        if cif is None:
            missing_struct += 1
            # still write zeros for these PDBs (so later scripts have complete coverage)
            rows = []
            for ch, fids in fasta_map[pdb].items():
                for fid in fids:
                    rows.append((fid, 0, 0, 0, 0))
            with SS_COUNTS_TSV.open("a", encoding="utf-8", newline="\n") as out:
                for fid, br, hr, bs, hx in rows:
                    out.write(f"{fid}\t{br}\t{hr}\t{bs}\t{hx}\n")
            with FAIL_LOG.open("a", encoding="utf-8") as f:
                f.write(f"[MISSING_STRUCT] {pdb}: no structure CIF filename contains this PDBID\n")
            continue

        jobs.append((pdb, cif, fasta_map[pdb]))

    print("===== DSSP SS COUNTS (robust) =====")
    print("Total PDBs in FASTA:", len(pdb_ids))
    print("Jobs with structure CIF found:", len(jobs))
    print("Missing structure CIF:", missing_struct)
    print("Workers:", WORKERS)
    print("Output:", SS_COUNTS_TSV.name)
    print("Fail log:", FAIL_LOG.name)

    ok = fail = 0
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process_one_pdb, pdb, cif, chains): (pdb, cif) for pdb, cif, chains in jobs}
        done = 0

        for fut in as_completed(futs):
            done += 1
            pdb, success, mode, rows, err = fut.result()

            if not success:
                fail += 1
                with FAIL_LOG.open("a", encoding="utf-8") as f:
                    f.write(f"[FAIL] {pdb} mode={mode}\n{err}\n{'-'*80}\n")
            else:
                ok += 1

            # append rows
            with tsv_lock:
                with SS_COUNTS_TSV.open("a", encoding="utf-8", newline="\n") as out:
                    for fid, br, hr, bs, hx in sorted(rows):
                        out.write(f"{fid}\t{br}\t{hr}\t{bs}\t{hx}\n")

            if done % 20 == 0 or done == len(jobs):
                elapsed = time.time() - start_time
                print(f"Progress: {done}/{len(jobs)} ok={ok} fail={fail} elapsed={elapsed:.1f}s")

    print("\nDONE")
    print("Wrote:", SS_COUNTS_TSV)
    print("Fail log:", FAIL_LOG)


if __name__ == "__main__":
    main()
# .\venv\Scripts\python.exe .\step_ss_counts_dssp_all.py
