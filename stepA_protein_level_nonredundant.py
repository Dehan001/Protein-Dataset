from __future__ import annotations
from pathlib import Path
from Bio import SeqIO

BASE = Path(__file__).resolve().parent
IN_FASTA = BASE / "outputs" / "all_chains.fasta"
OUT_FASTA = BASE / "outputs" / "protein_level_nonredundant.fasta"
OUT_LIST = BASE / "outputs" / "protein_level_selected_chains.txt"

def pdb_id_from_fasta_id(fid: str) -> str:
    # Example: 6K4M_EMD-9915_A  -> 6K4M
    return fid.split("_", 1)[0].upper()

def main():
    if not IN_FASTA.exists():
        raise FileNotFoundError(IN_FASTA)

    best = {}  # pdb_id -> (len, record)
    for rec in SeqIO.parse(str(IN_FASTA), "fasta"):
        fid = rec.id
        pdb = pdb_id_from_fasta_id(fid)
        L = len(rec.seq)
        if pdb not in best or L > best[pdb][0]:
            best[pdb] = (L, rec)

    selected = [v[1] for v in best.values()]
    selected.sort(key=lambda r: pdb_id_from_fasta_id(r.id))

    OUT_FASTA.parent.mkdir(exist_ok=True)
    SeqIO.write(selected, str(OUT_FASTA), "fasta")

    with OUT_LIST.open("w", encoding="utf-8", newline="\n") as f:
        f.write("pdb_id\tchosen_chain_fasta_id\tlength\n")
        for pdb, (L, rec) in sorted(best.items()):
            f.write(f"{pdb}\t{rec.id}\t{L}\n")

    print("Protein-level nonredundant chains:", len(selected))
    print("Wrote:", OUT_FASTA)
    print("Wrote:", OUT_LIST)

if __name__ == "__main__":
    main()
# .\venv\Scripts\python.exe .\stepA_protein_level_nonredundant.py
