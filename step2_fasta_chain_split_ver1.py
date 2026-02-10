import os
import re
import requests
from pathlib import Path

BASE = Path(__file__).resolve().parent
PDB_IDS_FILE = BASE / "pdb_ids.txt"
EM_MAPS_DIR = BASE / "em_maps"
FASTA_RAW_DIR = BASE / "fasta_raw"
CHAIN_FASTA_DIR = BASE / "chain_fasta"
OUT_DIR = BASE / "outputs"

FASTA_RAW_DIR.mkdir(exist_ok=True)
CHAIN_FASTA_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)

def emdb_ids_for_pdb(pdb_id: str):
    """
    Read EMDB IDs from downloaded map filenames:
    em_maps/<PDBID>_EMD-xxxxx.mrc  -> returns ["EMD-xxxxx", ...]
    """
    emdbs = []
    prefix = f"{pdb_id}_EMD-"
    for f in EM_MAPS_DIR.glob(f"{pdb_id}_EMD-*.mrc"):
        name = f.stem  # e.g., 6K4M_EMD-12345
        parts = name.split("_", 1)
        if len(parts) == 2 and parts[1].startswith("EMD-"):
            emdbs.append(parts[1])
    # Stable order
    return sorted(set(emdbs))

def download_fasta(pdb_id: str) -> str:
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text

def parse_chain_list_from_header(header: str):
    """
    Accept both '|Chain ...|' and '|Chains ...|'
    Extract 1-letter chain IDs even if written like 'A[auth 1]' or 'A [auth 1]'
    Keep ordering; caller will keep chains[0] (first 1-letter chain).
    """
    # Match either "|Chain ...|" or "|Chains ...|"
    m = re.search(r"\|Chains?\s+([^|]+)\|", header)
    if not m:
        return []

    chains_part = m.group(1)

    # Split by comma
    tokens = [t.strip() for t in chains_part.split(",")]

    one_letter = []
    for t in tokens:
        # Grab first letter before any bracket/space (e.g., "A[auth 1]" -> "A")
        m2 = re.match(r"^([A-Za-z])(?:\s*\[.*\])?$", t)
        if m2:
            one_letter.append(m2.group(1).upper())
            continue

        # Also handle cases like "A[auth 1]" but with extra text
        m3 = re.match(r"^([A-Za-z])\b", t)
        if m3:
            one_letter.append(m3.group(1).upper())

    # Keep only A–Z
    one_letter = [c for c in one_letter if re.fullmatch(r"[A-Z]", c)]
    return one_letter

    """
    Header example:
    >3C91_2|Chains AA[auth 1], BA[auth 2], H, I, J, K, ...|Protein name|Organism

    We need chain IDs. Rule from instructions:
    - If multiple identical chains merged into one FASTA entry, KEEP ONLY the FIRST chain having 1-letter name.
    - In practice, we treat 1-letter chains (A-Z) as "good" chain IDs.
    """
    # Find "...|Chains ...|"
    m = re.search(r"\|Chains\s+([^|]+)\|", header)
    if not m:
        return []

    chains_part = m.group(1)
    # Split by commas, strip whitespace
    tokens = [t.strip() for t in chains_part.split(",")]

    # Extract 1-letter chain IDs (A-Z)
    one_letter = []
    for t in tokens:
        # token might be "H" or "AA[auth 1]" etc.
        # We only accept exactly one uppercase letter
        if re.fullmatch(r"[A-Z]", t):
            one_letter.append(t)

    return one_letter  # could be [] if none

def split_fasta_records(fasta_text: str):
    """
    Returns list of (header_line_without_>, sequence_string)
    """
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

# Outputs:
# 1) chain_fasta/<PDBID>_<EMDBID>_<chain>.fasta (one per chain-level record we keep)
# 2) outputs/all_chains.fasta (combined, ready for clustering)
combined_path = OUT_DIR / "all_chains.fasta"
combined = open(combined_path, "w", newline="\n")

pdb_ids = [line.strip() for line in PDB_IDS_FILE.read_text().splitlines() if line.strip()]

kept = 0
no_emdb = 0
no_chain = 0

for pdb in pdb_ids:
    emdbs = emdb_ids_for_pdb(pdb)
    if not emdbs:
        # Your Step 1 filter should prevent this, but keep safe
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
        chains = parse_chain_list_from_header(header)

        # Instruction rule: if merged, keep only FIRST 1-letter chain
        if not chains:
            no_chain += 1
            continue

        chain_id = chains[0]  # keep first 1-letter chain only

        for emdb in emdbs:
            out_id = f"{pdb}_{emdb}_{chain_id}"
            out_file = CHAIN_FASTA_DIR / f"{out_id}.fasta"

            out_file.write_text(f">{out_id}\n{seq}\n", encoding="utf-8", newline="\n")
            combined.write(f">{out_id}\n{seq}\n")
            kept += 1

combined.close()

print("===== STEP 2 SUMMARY =====")
print(f"PDB entries processed: {len(pdb_ids)}")
print(f"Chain-level FASTA written: {kept}")
print(f"PDB with no EMDB (should be 0): {no_emdb}")
print(f"FASTA records with no 1-letter chain parsed: {no_chain}")
print(f"Combined FASTA for clustering: {combined_path}")
