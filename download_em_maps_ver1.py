import os, gzip, shutil
import requests

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_IDS_FILE = os.path.join(BASE_DIR, "pdb_ids.txt")
OUT_DIR = os.path.join(BASE_DIR, "em_maps")

os.makedirs(OUT_DIR, exist_ok=True)

def get_emdb_ids(pdb_id: str):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return []
    data = r.json()
    return data.get("rcsb_entry_container_identifiers", {}).get("emdb_ids", []) or []

def emdb_map_url(emdb: str):
    # emdb like "EMD-12345"
    num = emdb.replace("EMD-", "")
    return f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{num}/map/emd_{num}.map.gz"

def download_file(url: str, out_path: str):
    r = requests.get(url, stream=True, timeout=120)
    if r.status_code != 200:
        return False
    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)
    return True

def gunzip_to_mrc(gz_path: str, mrc_path: str):
    with gzip.open(gz_path, "rb") as f_in, open(mrc_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

with open(PDB_IDS_FILE, "r") as f:
    pdb_ids = [line.strip() for line in f if line.strip()]

downloaded = 0
skipped = 0
no_emdb = 0
failed = 0

for pdb in pdb_ids:
    emdbs = get_emdb_ids(pdb)
    if not emdbs:
        print(f"[NO EMDB] {pdb}")
        no_emdb += 1
        continue

    for emdb in emdbs:
        out_mrc = os.path.join(OUT_DIR, f"{pdb}_{emdb}.mrc")   # <PDBID>_<EMDBID>.mrc
        if os.path.exists(out_mrc):
            print(f"[SKIP] {os.path.basename(out_mrc)}")
            skipped += 1
            continue

        url = emdb_map_url(emdb)
        tmp_gz = out_mrc + ".map.gz"

        print(f"[DL] {pdb} -> {emdb}")
        ok = download_file(url, tmp_gz)
        if not ok:
            print(f"[FAIL] {pdb} {emdb} (map not found / server error)")
            failed += 1
            continue

        gunzip_to_mrc(tmp_gz, out_mrc)
        os.remove(tmp_gz)
        downloaded += 1

print("\n===== SUMMARY =====")
print(f"Downloaded maps: {downloaded}")
print(f"Skipped existing: {skipped}")
print(f"PDBs with no EMDB: {no_emdb}")
print(f"Failed downloads: {failed}")
