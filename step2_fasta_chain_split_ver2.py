import os
import re
import requests
from pathlib import Path

# ===================== PATHS =====================
BASE = Path(__file__).resolve().parent
PDB_IDS_FILE = BASE / "pdb_ids.txt"
EM_MAPS_DIR = BASE / "em_maps"
FASTA_RAW_DIR = BASE / "fasta_raw"
CHAIN_FASTA_DIR = BASE / "chain_fasta"
OUT_DIR = BASE / "outputs"

FASTA_RAW_DIR.mkdir(exist_ok=True)
CHAIN_FASTA_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)

# ===================== FILTER SETTINGS =====================
MAX_X_FRAC = 0.20  # Drop sequence if X > 20%

# ===================== HELPERS =====================
def x_fraction(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 1.0
    return seq.count("X") / len(seq)

# ===================== EMDB ↔ PDB =====================
def emdb_ids_for_pdb(pdb_id: str):
    """
    Read EMDB IDs from filenames:
    em_maps/<PDBID>_EMD-xxxxx.mrc
    """
    emdbs = []
    for f in EM_MAPS_DIR.glob(f"{pdb_id}_EMD-*.mrc"):
        name = f.stem
        parts = name.split("_", 1)
        if len(parts) == 2 and parts[1].startswith("EMD-"):
            emdbs.append(parts[1])
    return sorted(set(emdbs))

# ===================== FASTA =====================
def download_fasta(pdb_id: str) -> str:
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text

def split_fasta_records(fasta_text: str):
    records = []
    header = None
    seq_lines = []

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines)))
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)

    if header is not None:
        records.append((header, "".join(seq_lines)))

    return records

# ===================== CHAIN PARSING =====================
def parse_chain_list_from_header(header: str):
    """
    Accepts:
    |Chain A|
    |Chains A[auth 1], B, C|
    Returns ordered 1-letter chain IDs
    """
    m = re.search(r"\|Chains?\s+([^|]+)\|", header)
    if not m:
        return []

    tokens = [t.strip() for t in m.group(1).split(",")]

    chains = []
    for t in tokens:
        m1 = re.match(r"^([A-Za-z])(?:\s*\[.*\])?$", t)
        if m1:
            chains.append(m1.group(1).upper())
            continue
        m2 = re.match(r"^([A-Za-z])\b", t)
        if m2:
            chains.append(m2.group(1).upper())

    return [c for c in chains if re.fullmatch(r"[A-Z]", c)]

# ===================== OUTPUT =====================
combined_path = OUT_DIR / "all_chains.fasta"
combined = open(combined_path, "w", newline="\n")

pdb_ids = [l.strip() for l in PDB_IDS_FILE.read_text().splitlines() if l.strip()]

kept = 0
no_emdb = 0
no_chain = 0
dropped_high_x = 0

# ===================== MAIN LOOP =====================
for pdb in pdb_ids:
    emdbs = emdb_ids_for_pdb(pdb)
    if not emdbs:
        emdbs = ["EMD-UNKNOWN"]
        no_emdb += 1

    raw_path = FASTA_RAW_DIR / f"{pdb}.fasta"
    if not raw_path.exists():
        fasta_text = download_fasta(pdb)
        raw_path.write_text(fasta_text, encoding="utf-8", newline="\n")
    else:
        fasta_text = raw_path.read_text(encoding="utf-8")

    records = split_fasta_records(fasta_text)

    for header, seq in records:
        seq_u = seq.upper()

        # ---------- DROP BAD SEQUENCES ----------
        if x_fraction(seq_u) > MAX_X_FRAC:
            dropped_high_x += 1
            continue

        chains = parse_chain_list_from_header(header)
        if not chains:
            no_chain += 1
            continue

        chain_id = chains[0]  # rule: keep FIRST valid 1-letter chain

        for emdb in emdbs:
            out_id = f"{pdb}_{emdb}_{chain_id}"
            out_file = CHAIN_FASTA_DIR / f"{out_id}.fasta"

            out_file.write_text(
                f">{out_id}\n{seq_u}\n",
                encoding="utf-8",
                newline="\n"
            )

            combined.write(f">{out_id}\n{seq_u}\n")
            kept += 1

combined.close()

# ===================== SUMMARY =====================
print("===== STEP 2 SUMMARY =====")
print(f"PDB entries processed: {len(pdb_ids)}")
print(f"Chain-level FASTA written: {kept}")
print(f"Dropped (X > 20%): {dropped_high_x}")
print(f"PDB with no EMDB: {no_emdb}")
print(f"FASTA records with no chain parsed: {no_chain}")
print(f"Combined FASTA: {combined_path}")
