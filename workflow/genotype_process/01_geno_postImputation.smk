#!/usr/bin/env python3

import pandas as pd
import os

def OUT(*parts):
    return os.path.join(config["paths"]["genotype_outdir"], *[p for p in parts if p])

mergeName = config["final_vcf_name"]
POST_IMPUTATION_DIR = config["post_imputation_dir"]

onsuccess:
    print("Imputed VCFs processed successfully!")
    print(f"Final output: {OUT('vcf', mergeName + '_ALL_qc.vcf.gz')}")

rule all:
    input: 
        OUT('vcf', mergeName + '_ALL_qc.vcf.gz'),
        OUT('reports', mergeName + '_ALL_qc_VCFreport.txt')

rule unzip_files:
    output:
        expand(POST_IMPUTATION_DIR + "/chr{chr}.dose.vcf.gz", chr=range(1, 23))
    params:
        post_dir=POST_IMPUTATION_DIR,
        zip_password=config['zip_password']
    shell:
        """
        # Unzip all chromosome files with password
        for chr in {{1..22}}; do
            if [[ -f "{params.post_dir}/chr_${{chr}}.zip" ]]; then
                echo "Extracting chr_${{chr}}.zip"
                unzip -o -P "{params.zip_password}" "{params.post_dir}/chr_${{chr}}.zip" -d "{params.post_dir}/"
            fi
        done
        """

rule concat_vcf:
    input:
        rules.unzip_files.output
    output:
        OUT('vcf', mergeName + '_ALL.vcf.gz')
    params:
        post_dir=POST_IMPUTATION_DIR,
        temp_vcf=OUT('vcf', mergeName + '_ALL.vcf'),
        samtools_ver=config['samtools'],
        vcf_dir=OUT('vcf')
    log:
        err=OUT("logs", "concat_vcf.err"),
        out=OUT("logs", "concat_vcf.out")
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p {params.vcf_dir}

        # Create file list for concatenation
        vcf_list=""
        for chr in {{1..22}}; do
            vcf_list="$vcf_list {params.post_dir}/chr${{chr}}.dose.vcf.gz"
        done

        bcftools concat $vcf_list | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -o {params.temp_vcf}
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
