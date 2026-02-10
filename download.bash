mkdir pdb_models
cd pdb_models

# For linux
# while read pdb; do
#   wget https://files.rcsb.org/download/${pdb}.cif
# done < ../pdb_ids.txt

Get-Content ..\pdb_ids.txt | ForEach-Object {
    Invoke-WebRequest "https://files.rcsb.org/download/$_.cif" -OutFile "$_.cif"
}

cd path\to\Protein

mkdir em_maps

# Create a file: Protein\download_em_maps.py

mkdir fasta_raw
mkdir chain_fasta
mkdir outputs

# cd C:\Protein
py -m venv venv
# Created an block
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned

.\venv\Scripts\Activate.ps1

pip install --upgrade pip
Chunks done: 80 | passing_edges=870 | elapsed=0.25h
Chunks done: 160 | passing_edges=1411 | elapsed=0.52h
Chunks done: 240 | passing_edges=2087 | elapsed=0.71h
Chunks done: 320 | passing_edges=3435 | elapsed=1.09h
Chunks done: 400 | passing_edges=5608 | elapsed=1.36h
Chunks done: 480 | passing_edges=10069 | elapsed=1.58h
Chunks done: 560 | passing_edges=11161 | elapsed=1.93h

===== STEP 1–2 DONE (FULL) =====
Passing edges written: 14947
Wrote: C:\Users\fd02629\Desktop\Protein Experiment\outputs\similar_pairs.tsv
Elapsed: 2.25 hours


cd "C:\Users\fd02629\Desktop\Protein Experiment"

mkdir mmseqs_work -ErrorAction SilentlyContinue
mkdir mmseqs_tmp  -ErrorAction SilentlyContinue

$threads = 8  # or 10

mmseqs createdb .\outputs\all_chains.fasta .\mmseqs_work\chainsDB

conda --version

conda create -n mmseqs -c conda-forge -c bioconda mmseqs2 -y
conda activate mmseqs

# no mmseqs package found for windows
# Now trying for linux for windows

wsl --install -d Ubuntu

sudo apt update
sudo apt install -y mmseqs2
mmseqs --version

$MMSEQS = "C:\Users\fd02629\Tools\mmseqs"
$env:Path += ";$MMSEQS;$MMSEQS\bin"
mmseqs version

# Records

cd "C:\Users\fd02629\Desktop\Protein Experiment"

(Get-Content .\pdb_ids.txt | Sort-Object -Unique) | Set-Content .\outputs\pdb_ids_unique.txt
(Get-Content .\outputs\pdb_ids_unique.txt).Count


===== STEP 5.1 (DSSP -> beta STRANDS per FASTA chain) =====
Total CIFs: 1099
Workers: 8
Progress: 10/1099 (ok=10, fail=0, zero=2) elapsed=15.5s
Progress: 20/1099 (ok=20, fail=0, zero=3) elapsed=98.0s
Progress: 30/1099 (ok=30, fail=0, zero=4) elapsed=125.0s
Progress: 40/1099 (ok=40, fail=0, zero=5) elapsed=164.1s
Progress: 50/1099 (ok=50, fail=0, zero=8) elapsed=253.2s
Progress: 60/1099 (ok=60, fail=0, zero=9) elapsed=274.7s
Progress: 70/1099 (ok=70, fail=0, zero=12) elapsed=318.2s
Progress: 80/1099 (ok=80, fail=0, zero=14) elapsed=335.7s
Progress: 90/1099 (ok=90, fail=0, zero=15) elapsed=352.8s
Progress: 100/1099 (ok=100, fail=0, zero=19) elapsed=398.9s
Progress: 110/1099 (ok=110, fail=0, zero=22) elapsed=430.0s
Progress: 120/1099 (ok=120, fail=0, zero=23) elapsed=455.0s
Progress: 130/1099 (ok=130, fail=0, zero=26) elapsed=482.3s
Progress: 140/1099 (ok=140, fail=0, zero=27) elapsed=504.0s
Progress: 150/1099 (ok=150, fail=0, zero=30) elapsed=522.2s
Progress: 160/1099 (ok=159, fail=1, zero=33) elapsed=1057.2s
Progress: 170/1099 (ok=169, fail=1, zero=36) elapsed=1077.2s
Progress: 180/1099 (ok=179, fail=1, zero=40) elapsed=1097.2s
Progress: 190/1099 (ok=189, fail=1, zero=42) elapsed=1109.1s
Progress: 200/1099 (ok=199, fail=1, zero=45) elapsed=1143.6s
Progress: 210/1099 (ok=209, fail=1, zero=46) elapsed=1180.5s
Progress: 220/1099 (ok=219, fail=1, zero=48) elapsed=1231.0s
Progress: 230/1099 (ok=229, fail=1, zero=49) elapsed=1304.2s
Progress: 240/1099 (ok=239, fail=1, zero=51) elapsed=1467.1s
Progress: 250/1099 (ok=249, fail=1, zero=52) elapsed=1595.6s
Progress: 260/1099 (ok=259, fail=1, zero=54) elapsed=1620.8s
Progress: 270/1099 (ok=269, fail=1, zero=54) elapsed=1636.4s
Progress: 280/1099 (ok=279, fail=1, zero=55) elapsed=1663.7s
Progress: 290/1099 (ok=287, fail=3, zero=62) elapsed=1841.3s
Progress: 300/1099 (ok=297, fail=3, zero=64) elapsed=1904.6s
Progress: 310/1099 (ok=307, fail=3, zero=73) elapsed=1916.6s
Progress: 320/1099 (ok=317, fail=3, zero=83) elapsed=1923.9s
Progress: 330/1099 (ok=327, fail=3, zero=88) elapsed=1940.1s
Progress: 340/1099 (ok=337, fail=3, zero=89) elapsed=1963.7s
Progress: 350/1099 (ok=347, fail=3, zero=90) elapsed=1986.4s
Progress: 360/1099 (ok=357, fail=3, zero=93) elapsed=2069.8s
Progress: 370/1099 (ok=367, fail=3, zero=94) elapsed=2126.6s
Progress: 380/1099 (ok=377, fail=3, zero=96) elapsed=2216.4s
Progress: 390/1099 (ok=387, fail=3, zero=98) elapsed=2224.5s
Progress: 400/1099 (ok=395, fail=5, zero=104) elapsed=2571.8s
Progress: 410/1099 (ok=404, fail=6, zero=105) elapsed=2630.2s
Progress: 420/1099 (ok=412, fail=8, zero=107) elapsed=2644.3s
Progress: 430/1099 (ok=422, fail=8, zero=108) elapsed=2668.9s
Progress: 440/1099 (ok=432, fail=8, zero=109) elapsed=2680.0s
Progress: 450/1099 (ok=442, fail=8, zero=109) elapsed=2687.7s
Progress: 460/1099 (ok=452, fail=8, zero=109) elapsed=2696.7s
Progress: 470/1099 (ok=462, fail=8, zero=111) elapsed=2704.8s
Progress: 480/1099 (ok=472, fail=8, zero=113) elapsed=2714.0s
Progress: 490/1099 (ok=482, fail=8, zero=113) elapsed=2782.9s
Progress: 500/1099 (ok=492, fail=8, zero=114) elapsed=2835.2s
Progress: 510/1099 (ok=501, fail=9, zero=116) elapsed=2932.5s
Progress: 520/1099 (ok=511, fail=9, zero=119) elapsed=2953.8s
Progress: 530/1099 (ok=521, fail=9, zero=121) elapsed=2967.7s
Progress: 540/1099 (ok=531, fail=9, zero=122) elapsed=2994.4s
Progress: 550/1099 (ok=541, fail=9, zero=123) elapsed=3097.3s
Progress: 560/1099 (ok=551, fail=9, zero=123) elapsed=3142.8s
Progress: 570/1099 (ok=560, fail=10, zero=128) elapsed=3602.9s
Progress: 580/1099 (ok=570, fail=10, zero=130) elapsed=3640.8s
Progress: 590/1099 (ok=580, fail=10, zero=131) elapsed=3659.1s
Progress: 600/1099 (ok=590, fail=10, zero=132) elapsed=3704.5s
Progress: 610/1099 (ok=599, fail=11, zero=133) elapsed=3722.9s
Progress: 620/1099 (ok=609, fail=11, zero=133) elapsed=3748.7s
Progress: 630/1099 (ok=619, fail=11, zero=135) elapsed=3763.9s
Progress: 640/1099 (ok=629, fail=11, zero=137) elapsed=3780.5s
Progress: 650/1099 (ok=638, fail=12, zero=140) elapsed=3905.5s
Progress: 660/1099 (ok=648, fail=12, zero=140) elapsed=3971.2s
Progress: 670/1099 (ok=658, fail=12, zero=141) elapsed=4000.2s
Progress: 680/1099 (ok=668, fail=12, zero=142) elapsed=4009.4s
Progress: 690/1099 (ok=678, fail=12, zero=143) elapsed=4021.2s
Progress: 700/1099 (ok=688, fail=12, zero=145) elapsed=4040.5s
Progress: 710/1099 (ok=698, fail=12, zero=146) elapsed=4065.0s
Progress: 720/1099 (ok=708, fail=12, zero=149) elapsed=4073.3s
Progress: 730/1099 (ok=718, fail=12, zero=156) elapsed=4154.5s
Progress: 740/1099 (ok=728, fail=12, zero=162) elapsed=4165.3s
Progress: 750/1099 (ok=738, fail=12, zero=169) elapsed=4172.5s
Progress: 760/1099 (ok=747, fail=13, zero=171) elapsed=4559.2s
Progress: 770/1099 (ok=756, fail=14, zero=174) elapsed=4607.5s
Progress: 780/1099 (ok=766, fail=14, zero=174) elapsed=4681.4s
Progress: 790/1099 (ok=775, fail=15, zero=176) elapsed=4738.7s
Progress: 800/1099 (ok=785, fail=15, zero=179) elapsed=4859.8s
Progress: 810/1099 (ok=794, fail=16, zero=181) elapsed=4929.5s
Progress: 820/1099 (ok=803, fail=17, zero=182) elapsed=4994.0s
Progress: 830/1099 (ok=813, fail=17, zero=186) elapsed=5030.3s
Progress: 840/1099 (ok=822, fail=18, zero=190) elapsed=5189.9s
Progress: 850/1099 (ok=832, fail=18, zero=193) elapsed=5247.4s
Progress: 860/1099 (ok=842, fail=18, zero=194) elapsed=5516.9s
Progress: 870/1099 (ok=851, fail=19, zero=196) elapsed=5654.1s
Progress: 880/1099 (ok=861, fail=19, zero=198) elapsed=5974.8s
Progress: 890/1099 (ok=871, fail=19, zero=200) elapsed=6034.6s
Progress: 900/1099 (ok=881, fail=19, zero=201) elapsed=6079.6s
Progress: 910/1099 (ok=890, fail=20, zero=204) elapsed=6188.2s
Progress: 920/1099 (ok=900, fail=20, zero=207) elapsed=6200.5s
Progress: 930/1099 (ok=910, fail=20, zero=207) elapsed=6240.4s
Progress: 940/1099 (ok=920, fail=20, zero=210) elapsed=6271.7s
Progress: 950/1099 (ok=930, fail=20, zero=212) elapsed=6356.9s
Progress: 960/1099 (ok=940, fail=20, zero=213) elapsed=6378.3s
Progress: 970/1099 (ok=950, fail=20, zero=215) elapsed=6472.0s
Progress: 980/1099 (ok=960, fail=20, zero=218) elapsed=6617.4s
Progress: 990/1099 (ok=970, fail=20, zero=221) elapsed=7179.5s
Progress: 1000/1099 (ok=979, fail=21, zero=228) elapsed=7537.0s
Progress: 1010/1099 (ok=989, fail=21, zero=233) elapsed=7554.1s
Progress: 1020/1099 (ok=999, fail=21, zero=236) elapsed=7588.6s
Progress: 1030/1099 (ok=1009, fail=21, zero=236) elapsed=7617.2s
Progress: 1040/1099 (ok=1019, fail=21, zero=238) elapsed=7634.7s
Progress: 1050/1099 (ok=1029, fail=21, zero=239) elapsed=7645.2s
Progress: 1060/1099 (ok=1039, fail=21, zero=242) elapsed=7673.9s
Progress: 1070/1099 (ok=1049, fail=21, zero=245) elapsed=7702.8s
Progress: 1080/1099 (ok=1059, fail=21, zero=248) elapsed=7737.7s
Progress: 1090/1099 (ok=1069, fail=21, zero=252) elapsed=7771.4s
Progress: 1099/1099 (ok=1077, fail=22, zero=256) elapsed=8424.2s

===== STEP 5.1 DONE =====
OK: 1077 FAIL: 22 ZERO(all chains 0): 256
Wrote: C:\Users\fd02629\Desktop\Protein Experiment\outputs\beta_counts.tsv
Fail log: C:\Users\fd02629\Desktop\Protein Experiment\outputs\dssp_failures.log
Zero log: C:\Users\fd02629\Desktop\Protein Experiment\outputs\dssp_zero_counts.log

cd "C:\Users\fd02629\Desktop\Protein Experiment"
dir .\dssp_out | sort LastWriteTime -Descending | select -First 15

DONE
Wrote: C:\Users\fd02629\Desktop\Protein Experiment\outputs\ss_counts_20260203_112915.tsv
Fail log: C:\Users\fd02629\Desktop\Protein Experiment\outputs\dssp_failures_20260203_112915.log

===== MMseqs2 CLUSTER RECORDS SUMMARY =====
FASTA file: all_chains_ss_nonzero.fasta
MMseqs TSV: mmseqs_clusters_ss_nonzero.tsv

Total chains in FASTA: 3780
TSV lines read: 3780
TSV skipped (IDs not in FASTA): 0
TSV self-pairs rep==mem: 1732

Chains missing from TSV (added as singleton clusters): 0
Chains covered by clustering (after adding missing): 3780 (should equal 3780)

Total clusters (unique reps): 1732
Singleton clusters: 1086
Non-singleton clusters: 646
Chains in singleton clusters: 1086
Chains in size>=2 clusters: 2694

Wrote files:
 - records_chain_ids_20260203_225432.txt
 - records_assignment_20260203_225432.tsv
 - records_clusters_20260203_225432.tsv
 - records_reps_20260203_225432.txt
 - records_singleton_clusters_20260203_225432.tsv
 - records_nonsingle_clusters_20260203_225432.tsv
 - records_singleton_chains_20260203_225432.txt
 - records_chains_in_size_ge2_20260203_225432.txt
 - records_missing_from_tsv_20260203_225432.txt
 - records_cluster_size_distribution_20260203_225432.tsv
