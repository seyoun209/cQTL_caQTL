#!/bin/bash
#SBATCH -p general
#SBATCH -t 1-00:00:00
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH -J ld_filter
#SBATCH -o /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_filter_%A_%a.out
#SBATCH -e /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/ld_filter_%A_%a.err
#SBATCH --array=1-44

module load r/4.5.0

# Condition and chromosome arrays
conds=("pbs" "fnf")
chrs=({1..22})


# Compute job indices
total_jobs=$(( ${#conds[@]} * ${#chrs[@]} ))
job_index=$((SLURM_ARRAY_TASK_ID - 1))

cond_index=$(( job_index / 22 ))
chr_index=$(( job_index % 22 ))

cond=${conds[$cond_index]}
chr=${chrs[$chr_index]}

echo "Launching condition=$cond chr=$chr"

Rscript 03_c_snp_filter.r $chr $cond

