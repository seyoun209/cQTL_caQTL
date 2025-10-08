#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=50g
#SBATCH -t 1-00:00:00
#SBATCH -J LD_calc_all
#SBATCH -o /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_chr%a.out
#SBATCH -e /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_chr%a.err
#SBATCH --array=1-22


module load plink/1.90b6.21

# --- Inputs ---
CHR=${SLURM_ARRAY_TASK_ID}
REF="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/ALL/chr_all_snps/1000G.GRCh38.20190312_chr${CHR}"
SNP_FILE="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/caQTL/data/sig_response_bychr/chr${CHR}_sig_response_snps.txt"

# --- Output dirs ---
OUT_DIR="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_ALL/ld0"
OUT_DIR_05="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_ALL/ld05"
mkdir -p "${OUT_DIR_05}"
mkdir -p "${OUT_DIR}"

# --- Parameters ---
WINDOW_KB=1000
THREADS=4


# --- Loop over SNPs in this chromosome’s list ---
while read -r snp; do
  [[ -z "$snp" ]] && continue
  
  echo "${snp}"

  # LD r2 >= 0 (full)
  plink --bfile "${REF}" \
        --ld-snp "${snp}" \
        --ld-window 2000000 \
        --ld-window-kb ${WINDOW_KB} \
        --ld-window-r2 0 \
        --r2 \
        --threads ${THREADS} \
        --out "${OUT_DIR}/${snp}"

  # LD r2 >= 0.5
  plink --bfile "${REF}" \
        --ld-snp "${snp}" \
        --ld-window 2000000 \
        --ld-window-kb ${WINDOW_KB} \
        --ld-window-r2 0.5 \
        --r2 \
        --threads ${THREADS} \
        --out "${OUT_DIR_05}/${snp}"

done < "${SNP_FILE}"
