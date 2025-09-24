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
        '00_ATAC_Preprocess_Run.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/00_atac_preprocess.smk --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
    '01_ATAC_verifybamID_Run.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/01_01_vcf_VerifyBamID.smk  --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
    '01_02_ATAC_callingPeak.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/01_02_callingPeak.smk --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
            'rna' | 'run_RNA_preprocessing')
            snakemake -j 1 --unlock -s workflow/RNA_preprocess.smk --configfile "config/RNAconfig.yaml" --profile config/profile_slurm
            ;;
    '02_ATAC_wasp_Run.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/02_wasp.smk --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
    '02_02_ATAC_wasp_callingPeak_Run.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/02_02_wasp_callingPeak.smk --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
    '01_genoPipe_Run.sbatch' | '01_geno')
            snakemake -j 1 --unlock -s workflow/genotype_process/00_geno_filter.smk --configfile config/Geno_config.yaml --profile config/profile_slurm_geno
            ;;
    '02_genoPipe_postImputation_Run.sbatch' | '02_genoPipe_postImputation_Run')
            snakemake -j 1 --unlock -s workflow/genotype_process/01_geno_postImputation.smk --configfile config/Geno_config.yaml --profile config/profile_slurm_geno
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."
