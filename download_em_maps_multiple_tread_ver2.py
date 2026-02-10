import os
import gzip
import shutil
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep

# -------------------------------
# Configuration
# -------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_IDS_FILE = os.path.join(BASE_DIR, "pdb_ids.txt")
OUT_DIR = os.path.join(BASE_DIR, "em_maps")
MAX_WORKERS = 10        # Number of concurrent downloads
RETRIES = 3             # Retry failed downloads
CHUNK_SIZE = 1024 * 1024  # 1 MB chunks

os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------------
# Helper functions
# -------------------------------

def get_emdb_ids(pdb_id: str):
    """Return list of EMDB IDs linked to a PDB ID."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            return []
        data = r.json()
        return data.get("rcsb_entry_container_identifiers", {}).get("emdb_ids", []) or []
    except Exception:
        return []

def emdb_map_url(emdb: str):
    """Return download URL for EMDB map."""
    num = emdb.replace("EMD-", "")
    return f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{num}/map/emd_{num}.map.gz"

def download_file(url: str, out_path: str):
    """Download a file with retries."""
    for attempt in range(1, RETRIES+1):
        try:
            r = requests.get(url, stream=True, timeout=120)
            if r.status_code != 200:
                raise Exception(f"HTTP {r.status_code}")
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                    if chunk:
                        f.write(chunk)
            return True
        except Exception as e:
            print(f"[WARN] Attempt {attempt} failed: {url} ({e})")
            sleep(2)  # brief pause before retry
    return False

def gunzip_to_mrc(gz_path: str, mrc_path: str):
    """Decompress .gz file to .mrc"""
    with gzip.open(gz_path, "rb") as f_in, open(mrc_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

def download_emdb_map(pdb_id: str):
    """Download all EMDB maps for a given PDB ID."""
    results = []
    emdbs = get_emdb_ids(pdb_id)
    if not emdbs:
        results.append((pdb_id, None, "no_emdb"))
        return results

    for emdb in emdbs:
        out_mrc = os.path.join(OUT_DIR, f"{pdb_id}_{emdb}.mrc")
        if os.path.exists(out_mrc):
            results.append((pdb_id, emdb, "skipped"))
            continue

        tmp_gz = out_mrc + ".map.gz"
        url = emdb_map_url(emdb)

        print(f"[DL] {pdb_id} -> {emdb}")
        if download_file(url, tmp_gz):
            gunzip_to_mrc(tmp_gz, out_mrc)
            os.remove(tmp_gz)
            results.append((pdb_id, emdb, "downloaded"))
        else:
            results.append((pdb_id, emdb, "failed"))
    return results

# -------------------------------
# Main
# -------------------------------

def main():
    with open(PDB_IDS_FILE, "r") as f:
        pdb_ids = [line.strip() for line in f if line.strip()]

    summary = {
        "downloaded": 0,
        "skipped": 0,
        "no_emdb": 0,
        "failed": 0
    }

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(download_emdb_map, pdb): pdb for pdb in pdb_ids}
        for future in as_completed(futures):
            for pdb, emdb, status in future.result():
                if status == "downloaded":
                    summary["downloaded"] += 1
                elif status == "skipped":
                    summary["skipped"] += 1
                elif status == "no_emdb":
                    summary["no_emdb"] += 1
                elif status == "failed":
                    summary["failed"] += 1
                print(f"[{status.upper()}] {pdb} {emdb if emdb else ''}")

    print("\n===== SUMMARY =====")
    print(f"Downloaded maps: {summary['downloaded']}")
    print(f"Skipped existing: {summary['skipped']}")
    print(f"PDBs with no EMDB: {summary['no_emdb']}")
    print(f"Failed downloads: {summary['failed']}")

if __name__ == "__main__":
    main()
