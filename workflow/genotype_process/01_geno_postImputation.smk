#!/usr/bin/env python3

import pandas as pd
import os

def OUT(*parts):
    return os.path.join(config["paths"]["genotype_outdir"], *[p for p in parts if p])

## Read in samplesheet
genos = pd.read_csv(config["geno"], sep=",")
genos = genos.astype(str)

## Read in samplesheet with imputed vcfs
vcfs = pd.read_csv(config["vcf"], sep=",")
vcfs = vcfs.astype(str)

## Grab vcf paths for all chromosomes
vcf_paths = vcfs['vcf_path'].tolist()

mergeName = '_'.join(genos.Proj.unique()) + '_' + '_'.join(genos['Batch'].unique())

onsuccess:
    print("Imputed VCFs processed successfully!")
    print(f"QC'd VCFs ready for allelic imbalance or QTL mapping pipelines.")
    print(f"Final output: {OUT('vcf', mergeName + '_ALL_qc.vcf.gz')}")

rule all:
    input: 
        OUT('vcf', mergeName + '_ALL_qc.vcf.gz'),
        OUT('reports', mergeName + '_ALL_qc_VCFreport.txt')

rule concat_vcf:
    input: 
        vcf_paths
    output:
        OUT('vcf', mergeName + '_ALL.vcf.gz')
    params:
        prefix=mergeName,
        vcf_dir=OUT('vcf'),
        temp_vcf=OUT('vcf', mergeName + '_ALL.vcf'),
        samtools_ver=config['samtools']
    log:
        err=OUT("logs", "concat_vcf.err"),
        out=OUT("logs", "concat_vcf.out")
    threads: 4
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p {params.vcf_dir}
        
        bcftools concat -o {params.temp_vcf} {input} 1> {log.out} 2> {log.err}
        bgzip {params.temp_vcf}
        """      

rule qc_vcf:
    input:
        OUT('vcf', mergeName + '_ALL.vcf.gz')
    output:
        OUT('vcf', mergeName + '_ALL_qc.vcf.gz')
    params:
        prefix=mergeName,
        geno=config['miss_call'],
        hwe=config['hwe'],
        maf=config['maf'],
        vcf_dir=OUT('vcf'),
        temp_prefix=OUT('vcf', mergeName + '_ALL_qc_temp'),
        temp_vcf=OUT('vcf', mergeName + '_ALL_qc.vcf'),
        filtered_variants='filtered_variants.txt',
        plink_ver=config['plink'],
        samtools_ver=config['samtools']
    threads: 4
    log:
        out=OUT("logs", mergeName + "_qc_vcf.out"),
        err=OUT("logs", mergeName + "_qc_vcf.err")
    shell:
        """
        module load plink/{params.plink_ver}
        module load samtools/{params.samtools_ver}
        mkdir -p {params.vcf_dir}
        
        # Run PLINK QC
        plink --vcf {input} --threads {threads} --const-fid 0 \
              --geno {params.geno} --hwe {params.hwe} --maf {params.maf} \
              --recode vcf --out {params.temp_prefix} 1> {log.out} 2> {log.err}
        
        # Move log file
        if [ -f {params.temp_prefix}.log ]; then
            mv {params.temp_prefix}.log {log.err}
        fi
        
        # Extract variant IDs that passed QC
        awk '!/^#/{{print $3}}' {params.temp_prefix}.vcf > {params.filtered_variants}
        
        # Filter original VCF to keep only QC-passed variants
        bcftools view -i 'ID=@{params.filtered_variants}' {input} \
                     -o {params.temp_vcf} -O v 2>> {log.err}
        
        # Clean up temporary files
        rm -f {params.temp_prefix}.nosex {params.filtered_variants}
        
        # Compress final VCF
        bgzip {params.temp_vcf}
        """

rule vcf_stats:
    input:
        OUT('vcf', mergeName + '_ALL_qc.vcf.gz')
    output:
        OUT('reports', mergeName + '_ALL_qc_VCFreport.txt')
    params:
        prefix=mergeName,
        reports_dir=OUT('reports'),
        samtools_ver=config['samtools']
    log:
        err=OUT("logs", "vcf_stats.err")
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p {params.reports_dir}
        
        bcftools stats {input} > {output} 2> {log.err}
        """
