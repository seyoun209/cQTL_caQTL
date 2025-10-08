#!/bin/bash
#SBATCH --job-name=ld_filter
#SBATCH -t 1-00:00:00
#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH --output=/work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_filter_%a.out
#SBATCH --error=/work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_filter_%a.err
#SBATCH --array=1-22

module load r/4.5.0

Rscript 00_ld_filter.r ${SLURM_ARRAY_TASK_ID}
