import os
import re
import requests
from pathlib import Path
import certifi

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
    r = requests.get(url, timeout=60, verify=certifi.where())
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
def parse_chain_tokens_from_header(header: str):
    """
    Extract chain tokens from FASTA header.

    Accepts headers like:
      |Chain A|
      |Chains AA[auth 1], BA[auth 2], H, I|
      |Chains A, B, C|

    Returns ordered chain IDs (strings) after removing [auth ...] and spaces.
    Examples: ["AA","BA","H","I"]
    """
    m = re.search(r"\|Chains?\s+([^|]+)\|", header)
    if not m:
        return []

    raw = m.group(1)

    tokens = []
    for t in raw.split(","):
        t = t.strip()
        # remove bracket annotation like [auth 1]
        t = re.sub(r"\s*\[.*\]\s*", "", t)
        # remove spaces
        t = t.replace(" ", "")
        if t:
            tokens.append(t.upper())

    return tokens

def choose_chain_id(chain_ids):
    """
    Always choose at least one chain from a duplicate-group list.
    Preference: shortest length first (1-letter preferred), then 2-letter, etc.
    Keeps original order within the same length.
    """
    if not chain_ids:
        return None

    # Keep only alphanumeric chain IDs (common in PDB/author chains)
    cleaned = [c for c in chain_ids if re.fullmatch(r"[A-Z0-9]+", c)]
    if not cleaned:
        return None

    best_len = min(len(c) for c in cleaned)
    for c in cleaned:
        if len(c) == best_len:
            return c
    return cleaned[0]

# ===================== OUTPUT =====================
combined_path = OUT_DIR / "all_chains.fasta"
combined = open(combined_path, "w", newline="\n", encoding="utf-8")

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
    try:
        if not raw_path.exists():
            fasta_text = download_fasta(pdb)
            raw_path.write_text(fasta_text, encoding="utf-8", newline="\n")
        else:
            fasta_text = raw_path.read_text(encoding="utf-8")
    except Exception as e:
        print(f"[WARN] FASTA download/read failed for {pdb}: {e}")
        continue

    records = split_fasta_records(fasta_text)

    for header, seq in records:
        seq_u = seq.upper()

        # ---------- DROP BAD SEQUENCES ----------
        if x_fraction(seq_u) > MAX_X_FRAC:
            dropped_high_x += 1
            continue

        # Parse chains from header (1-letter preferred, else 2-letter, etc.)
        chain_ids = parse_chain_tokens_from_header(header)
        chain_id = choose_chain_id(chain_ids)

        # If still nothing, skip (rare; header didn't contain chain list)
        if not chain_id:
            no_chain += 1
            continue

        # Write one chain-level sequence per EMDB ID
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
