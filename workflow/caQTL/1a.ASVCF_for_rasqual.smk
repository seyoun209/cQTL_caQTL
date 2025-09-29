#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

# Load and process samples to get MAIN_PBS and MAIN_FNF
samples = pd.read_csv(config["samplesheet"], sep='\t')
samples = samples.astype(str)
s_config = config['sample_logic']

# Handle replicate & femur naming logic
donor_col = s_config['femur_donor_col']
if str(s_config.get('Sample_handing', 'false')).lower() == 'true':
    rep_col = s_config['replicate_column']
    is_replicate = samples[rep_col].str.contains('replicate', case=False, na=False)
    
    tissue_col = s_config['femur_column']
    is_femur = samples[tissue_col].str.contains('Femur', na=False) | \
               samples[donor_col].str.endswith('FE', na=False)
    
    # Create full sample names
    full_df = samples.copy()
    full_df['replicate_suffix'] = np.where(is_replicate, 'replicate', '1')
    full_df[donor_col] = np.where(
        is_femur & ~full_df[donor_col].str.endswith('FE'),
        full_df[donor_col] + 'FE',
        full_df[donor_col]
    )
    full_name_cols = s_config['base_name_columns'] + ['replicate_suffix']
    samples['full_sample_name'] = full_df[full_name_cols].astype(str).agg('_'.join, axis=1)
    
    is_special = is_replicate | is_femur
else:
    std_name = samples[s_config['base_name_columns']].astype(str).agg('_'.join, axis=1) + '_1'
    samples['full_sample_name'] = std_name
    is_special = pd.Series(False, index=samples.index)

# Define sample groups
group_col = s_config['condition_column']
MAIN_SAMPLES = samples[~is_special]['full_sample_name'].tolist()
MAIN_PBS = samples[(samples[group_col] == s_config['control_group_value']) & ~is_special]['full_sample_name'].tolist()
MAIN_FNF = samples[(samples[group_col] == s_config['treatment_group_value']) & ~is_special]['full_sample_name'].tolist()

# Define the VCF groups  
VCF_GROUPS = ["pbs_geno", "fnf_geno"]
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

onsuccess:
    print("RASQUAL VCF preparation completed successfully!")

rule all:
    input:
        # 01. MAF filtered VCFs
        expand(OUT("geno", "{group}", "01_maf_filter", "finalFiltered.vcf.gz"), group=VCF_GROUPS),
        expand(OUT("geno", "{group}", "01_maf_filter", "finalFiltered.vcf.gz.tbi"), group=VCF_GROUPS),
        expand(OUT("geno", "{group}", "01_maf_filter", "finalFiltered_stats.txt"), group=VCF_GROUPS),
        # 02. PCA results
        expand(OUT("geno", "{group}", "02_pca", "cqtl.eigenval"), group=VCF_GROUPS),
        expand(OUT("geno", "{group}", "02_pca", "cqtl.eigenvec"), group=VCF_GROUPS),
        # 03. RASQUAL AS-VCFs by chromosome
        expand(OUT("geno", "{group}", "03_asvcf", "{chr}.{group}.ASCounts.vcf.gz"), 
               group=VCF_GROUPS, chr=CHROMOSOMES),
        # 04. eigenMT recode012 by chromosome
        expand(OUT("geno", "{group}", "04_recode", "{chr}.{group}.traw"),
               group=VCF_GROUPS, chr=CHROMOSOMES)

rule filter_by_maf_vcf:
    input:
        vcf=OUT("geno", "{group}", "{group}.vcf.gz"),
        tbi=OUT("geno", "{group}", "{group}.vcf.gz.tbi")
    output:
        filtered_vcf=OUT("geno", "{group}", "01_maf_filter", "finalFiltered.vcf.gz"),
        filtered_tbi=OUT("geno", "{group}", "01_maf_filter", "finalFiltered.vcf.gz.tbi")
    params:
        bcftools_ver=config['bcftools'],
        samtools_ver=config['samtoolsVers'],
        maf=config['maf_threshold']
    log:
        OUT("logs", "rasqual_maf_filter_{group}.log")
    shell:
        """
        module load bcftools/{params.bcftools_ver}
        module load samtools/{params.samtools_ver}
        mkdir -p $(dirname {output.filtered_vcf})
        
        # Filter by MAF - no sample name changes needed
        bcftools view {input.vcf} -i 'MAF >= {params.maf}' -Oz -o {output.filtered_vcf} 2> {log}
        
        # Index the final VCF
        tabix -p vcf {output.filtered_vcf} 2>> {log}
        """

rule generate_stats:
    input:
        vcf=rules.filter_by_maf_vcf.output.filtered_vcf
    output:
        stats=OUT("geno", "{group}", "01_maf_filter", "finalFiltered_stats.txt")
    params:
        bcftools_ver=config['bcftools']
    log:
        OUT("logs", "rasqual_stats_{group}.log")
    shell:
        """
        module load bcftools/{params.bcftools_ver}
        
        # Generate stats
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        """

rule create_pca:
    input:
        vcf=rules.filter_by_maf_vcf.output.filtered_vcf
    output:
        eigenval=OUT("geno", "{group}", "02_pca", "cqtl.eigenval"),
        eigenvec=OUT("geno", "{group}", "02_pca", "cqtl.eigenvec")
    params:
        plink_ver=config['plink'],
        output_prefix=lambda wc: OUT("geno", wc.group, "02_pca", "cqtl")
    log:
        OUT("logs", "rasqual_pca_{group}.log")
    shell:
        """
        module load plink/{params.plink_ver}
        mkdir -p $(dirname {output.eigenval})
        
        # Create PCA
        plink --vcf {input.vcf} \
              --pca \
              --out {params.output_prefix} \
              --double-id 2> {log}
        """

rule create_bam_lists:
    output:
        bam_list=OUT("geno", "{group}", "03_asvcf", "{group}_bamlist.txt")
    params:
        bam_dir=config['bam_dir']
    run:
        import os
        
        # Determine which samples to use based on group
        if wildcards.group == "pbs_geno":
            sample_list = MAIN_PBS
        else:  # fnf_geno
            sample_list = MAIN_FNF
            
        # Create output directory
        os.makedirs(os.path.dirname(output.bam_list), exist_ok=True)
        
        # Find corresponding BAM files
        bam_files = []
        for sample in sample_list:
            bam_path = os.path.join(params.bam_dir, f"{sample}.sorted_final.bam")
            if os.path.exists(bam_path):
                bam_files.append(bam_path)
        
        # Write BAM list
        with open(output.bam_list, 'w') as f:
            for bam_file in sorted(bam_files):
                f.write(bam_file + '\n')
        
        print(f"Created BAM list for {wildcards.group}: {len(bam_files)} files")

rule create_asvcf:
    input:
        vcf=rules.filter_by_maf_vcf.output.filtered_vcf,
        bam_list=rules.create_bam_lists.output.bam_list
    output:
        asvcf=OUT("geno", "{group}", "03_asvcf", "{chr}.{group}.ASCounts.vcf.gz"),
        asvcf_tbi=OUT("geno", "{group}", "03_asvcf", "{chr}.{group}.ASCounts.vcf.gz.tbi")
    params:
        samtools_ver=config['samtoolsVers'],
        rasqual_dir=config['rasqual_dir'],
        temp_dir=lambda wc: OUT("geno", wc.group, "03_asvcf", "temp")
    log:
        OUT("logs", "create_asvcf_{group}_{chr}.log")
    shell:
        """
        module load samtools/{params.samtools_ver}
        export RASQUALDIR={params.rasqual_dir}
        
        mkdir -p {params.temp_dir}
        mkdir -p $(dirname {output.asvcf})
        
        temp_vcf="{params.temp_dir}/temp_{wildcards.chr}.vcf.gz"
        
        # Extract chromosome
        bcftools view {input.vcf} {wildcards.chr} -Oz -o $temp_vcf 2> {log}
        
        # Index temp VCF
        tabix -p vcf $temp_vcf 2>> {log}
        
        # Run createASVCF
        $RASQUALDIR/src/ASVCF/createASVCF.sh paired_end {input.bam_list} $temp_vcf {output.asvcf} atac2 >> {log}
        
        # Index final AS-VCF
        tabix -p vcf {output.asvcf} 2>> {log}
        
        # Clean up temp files
        rm -f $temp_vcf $temp_vcf.tbi
        """

rule recode_vcf_for_eigenmt:
    input:
        vcf=rules.create_asvcf.output.asvcf
    output:
        traw=OUT("geno", "{group}", "04_recode", "{chr}.{group}.traw"),
    params:
        plink2_ver=config['plink2'],
        out_prefix=lambda wc: OUT("geno", wc.group, "04_recode", f"{wc.chr}.{wc.group}")
    threads: 4
    resources:
        mem_mb=8000,
        time="1:00:00"
    log:
        OUT("logs", "recode_{group}_{chr}.log")
    shell:
        """
        module load plink/{params.plink2_ver}
        mkdir -p $(dirname {output.traw})
        
        # Extract chromosome and recode to A-transpose format
        plink2 --vcf {input.vcf} \
               --chr {wildcards.chr} \
               --const-fid 0 \
               --recode A-transpose \
               --keep-allele-order \
               --out {params.out_prefix} 2> {log}
        
        # Fix header: remove "0_" prefix from sample IDs
        # The sed command removes "0_" after tab or at start of line
        sed -i '1s/\\(^\\|\\t\\)0_/\\1/g' {output.traw}
        """
