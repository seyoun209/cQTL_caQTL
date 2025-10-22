#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16g
#SBATCH -t 2:00:00
#SBATCH -J calculate_LD
#SBATCH -o /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_chr%a.out
#SBATCH -e /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_chr%a.err
#SBATCH --array=1-22

module load plink/1.90b6.21

# Make as the bfile for the vcf file
CHR=${SLURM_ARRAY_TASK_ID}
#VCF="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/pbs_geno/03_asvcf/chr${CHR}.pbs_geno.ASCounts.vcf.gz"
#OUT="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/pbs_geno/05_bfile_asvcf/chr${CHR}"
#mkdir -p /work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/pbs_geno/05_bfile_asvcf

#plink --vcf ${VCF} --const-fid --make-bed --out ${OUT}


BFILE="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/pbs_geno/05_bfile_asvcf/chr${CHR}"
SNP_LIST="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/caQTL/data/03_LD_snplist/window_10kb/chr${CHR}_SNPList.txt"
OUTPUT_DIR="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/caQTL/data/04_LD_results/window_10kb/ld0"
mkdir -p ${OUTPUT_DIR}


# Calculate LD for all
plink --bfile ${BFILE} \
      --extract ${SNP_LIST} \
      --ld-snp-list ${SNP_LIST} \
      --r2 \
      --ld-window-kb 1000 \
      --ld-window-r2 0 \
      --ld-window 2000000 \
      --out ${OUTPUT_DIR}/chr${CHR}.LD

echo "Chr ${CHR} LD calculation complete"
