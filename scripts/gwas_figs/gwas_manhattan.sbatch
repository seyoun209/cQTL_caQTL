#!/bin/bash
#SBATCH -p general
#SBATCH -t 1-00:00:00
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH -J gwas_fig
#SBATCH -o /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/gwas_fig_%A_%a.out
#SBATCH -e /work/users/s/e/seyoun/cQTL_caQTL/atac_output/logs/gwas_fig_%A_%a.err

module load r/4.5.0

Rscript gwas_manhattan.r
