#! /bin/bash -login

## Exit if any command fails
set -e

## Load required modules
module load python/3.12.4

## Activate virtual environment
source env/bin/activate

case $1 in

    '-h' | '--help' | ' ' | '')
            echo -e '\e[31mSpecify which workflow to unlock (i.e. run_RNAprocessing)'
            exit 2
            ;;
        '01_ATAC_Preprocess_Run.sbatch')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/ATAC_preatac_preprocess.smk --configfile "config/ATAC_config.yaml" --profile config/profile_slurm
            ;;
    'run_VCFpreprocessing')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/ATAC_preatac_preprocess.smk --configfile "config/ATAC_config.yaml" --profile config/profile_slurm
            ;;
    'run_SubsetPreprocessing')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/ATAC_preatac_preprocess.smk --configfile "config/ATAC_config.yaml" --profile config/profile_slurm
            ;;
            'rna' | 'run_RNA_preprocessing')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/RNA_preprocess.smk --configfile "config/RNAconfig.yaml" --profile config/profile_slurm
            ;;
    'wasp' | 'run_WASP_ATAC')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/WASP_atac.smk --configfile "config/ATACconfig.yaml" --profile config/profile_slurm
            ;;
    'signal' | 'run_signalTrack')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/signal_track.smk --configfile "config/ATACconfig.yaml" --profile config/profile_slurm
            ;;
    'peak' | 'run_peakcalling')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/peak_counting.smk --configfile "config/ATACconfig.yaml" --profile config/profile_slurm
            ;;
    'rna_signal' | 'run_RNAsignal')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s workflow/RNA_signal_track.smk --configfile "config/RNAconfig.yaml" --profile config/profile_slurm
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."
