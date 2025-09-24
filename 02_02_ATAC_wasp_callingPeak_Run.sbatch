#! /bin/bash -login
#SBATCH -J atac_wasp_peakcalling
#SBATCH -t 11-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p general
#SBATCH --mem=4gb
#SBATCH -o workflow_logs/%x-%j.out

## Exit if any command fails
set -e

## Load required modules
module load python/3.12.4

## Create and activate virtual environment with requirements
python3 -m venv env
source env/bin/activate && pip3 install -r config/requirements.txt
#pip install --upgrade pip

## Make directory for slurm logs
mkdir -p atac_output/logs_slurm workflow_logs


## Execute common preprocess
snakemake -s workflow/atac_process/02_02_wasp_callingPeak.smk \
  --configfile config/ATAC_config.yaml \
  --profile config/profile_slurm_atac \
  -j 50 --max-jobs-per-second 50 --max-status-checks-per-second 1 \
  --rerun-incomplete -p --latency-wait 200

echo "ATAC wasp peakcalling done!"
echo "Output files are in: atac_output/"


