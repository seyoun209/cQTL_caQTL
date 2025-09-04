#! /bin/bash -login
#SBATCH -J atac_preprocess
#SBATCH -t 11-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2  
#SBATCH -p general
#SBATCH --mem=4gb
#SBATCH -o workflow_logs/genoPipe_%x-%j.out

## Exit if any command fails
set -e

## Load required modules
module load python/3.12.4

## Create and activate virtual environment with requirements
python3 -m venv env
source env/bin/activate && pip3 install -r config/requirements.txt
#pip install --upgrade pip

## Make directory for slurm logs
mkdir -p geno_output/logs_slurm workflow_logs



