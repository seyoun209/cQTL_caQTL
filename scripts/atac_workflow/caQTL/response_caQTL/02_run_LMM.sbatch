#!/bin/bash
#SBATCH -p general
#SBATCH --mem=40g
#SBATCH -t 3-00:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -J response_qtl
#SBATCH -o /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/LMM_%a.out
#SBATCH -e /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/LMM_%a.err
#SBATCH --array=1-22

module load r/4.5.0

echo "Running LMM caQTL job for chromosome ${SLURM_ARRAY_TASK_ID}"

Rscript 02_run_LMM.r ${SLURM_ARRAY_TASK_ID}

echo "Finished chromosome ${SLURM_ARRAY_TASK_ID}"
