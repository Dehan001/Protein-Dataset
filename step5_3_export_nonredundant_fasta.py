from __future__ import annotations

from pathlib import Path
from typing import Set, List
from Bio import SeqIO  # pip install biopython

# ============== SETTINGS ==============
BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"

ALL_CHAINS_FASTA = OUT / "all_chains.fasta"
NONRED_IDS = OUT / "nonredundant_rep_names.txt"

OUT_FASTA = OUT / "nonredundant_reps.fasta"
MISSING_TXT = OUT / "nonredundant_missing_ids.txt"
# =====================================


def load_ids(path: Path) -> List[str]:
    ids: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                ids.append(s)
    return ids


def main():
    if not ALL_CHAINS_FASTA.exists():
        raise FileNotFoundError(f"Missing {ALL_CHAINS_FASTA}")
    if not NONRED_IDS.exists():
        raise FileNotFoundError(f"Missing {NONRED_IDS}. Run step5_2_select_representatives.py first.")

    wanted_list = load_ids(NONRED_IDS)
    wanted: Set[str] = set(wanted_list)

    written = 0
    found: Set[str] = set()

    # Stream-read FASTA and write only selected records
    with OUT_FASTA.open("w", encoding="utf-8", newline="\n") as out_f:
        for rec in SeqIO.parse(str(ALL_CHAINS_FASTA), "fasta"):
            rid = rec.id.strip()
            if rid in wanted:
                SeqIO.write(rec, out_f, "fasta")
                written += 1
                found.add(rid)

    missing = [x for x in wanted_list if x not in found]  # preserve original order
    MISSING_TXT.write_text("\n".join(missing) + ("\n" if missing else ""), encoding="utf-8")

    print("===== STEP 5.3 DONE =====")
    print("Wanted IDs:", len(wanted_list))
    print("Unique wanted IDs:", len(wanted))
    print("Written to FASTA:", written)
    print("Output FASTA:", OUT_FASTA)

    if missing:
        print("WARNING: Some IDs were not found in all_chains.fasta")
        print("Missing count:", len(missing))
        print("Missing list saved to:", MISSING_TXT)
        print("First 10 missing IDs:", missing[:10])
    else:
        print("All IDs were found. ✅")
        # Keep file empty to indicate no missing
        if MISSING_TXT.exists():
            # already written empty content above
            pass


if __name__ == "__main__":
    main()

# Run:
# .\venv\Scripts\python.exe -m pip install biopython
# .\venv\Scripts\python.exe .\step5_3_export_nonredundant_fasta.py
