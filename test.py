from pathlib import Path

p = Path("dssp_out/8GHA.dssp")
lines = p.read_text(errors="ignore").splitlines()

# find residue table
start = None
for i, ln in enumerate(lines):
    if ln.startswith("  #  RESIDUE"):
        start = i+1
        break

print("Found residue table:", start is not None)

if start:
    chains = {}
    for ln in lines[start:]:
        if len(ln) < 17: 
            continue
        ch = ln[11].strip()
        ss = ln[16].strip()
        if not ch: 
            continue
        chains.setdefault(ch, 0)
        if ss in {"E","B","H","G","I"}:
            chains[ch] += 1
    print("Chain counts (any E/B/H/G/I residues):", chains)
