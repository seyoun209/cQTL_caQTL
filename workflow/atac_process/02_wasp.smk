#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

#------------------------------------------------------
# Load and process samples
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
    
    # Create full sample names (used in BAMs)
    full_df = samples.copy()
    full_df['replicate_suffix'] = np.where(is_replicate, 'replicate', '1')
    full_df[donor_col] = np.where(
        is_femur & ~full_df[donor_col].str.endswith('FE'),
        full_df[donor_col] + 'FE',
        full_df[donor_col]
    )
    full_name_cols = s_config['base_name_columns'] + ['replicate_suffix']
    samples['full_sample_name'] = full_df[full_name_cols].astype(str).agg('_'.join, axis=1)
else:
    std_name = samples[s_config['base_name_columns']].astype(str).agg('_'.join, axis=1) + '_1'
    samples['full_sample_name'] = std_name

# Filter out excluded samples if needed
exclude_options = config.get('exclude_options', {})
if exclude_options.get('exclude_replicates', False):
    samples = samples[~is_replicate]
if exclude_options.get('exclude_synovium', False):
    samples = samples[~samples['Tissue'].str.contains('synovium', case=False, na=False)]
if exclude_options.get('exclude_femur', False):
    samples = samples[~is_femur]

# Group samples by condition
group_col = s_config['condition_column']
condition_samples = samples.groupby(group_col)['full_sample_name'].apply(list).to_dict()

# Create BAM file mappings
BAM_FILES = {}
BAI_FILES = {}
for _, row in samples.iterrows():
    sample_name = row['full_sample_name']
    bam_path = f"{config['paths']['base_dir']}/atac_output/filtered/blk_filter/{sample_name}.sorted_final.bam"
    BAM_FILES[sample_name] = bam_path
    BAI_FILES[sample_name] = bam_path + ".bai"

ALL_SAMPLES = list(BAM_FILES.keys())
VCF_GROUPS = ["pbs_geno", "fnf_geno", "pbs_special", "fnf_special"]
CHROMOSOMES = list(range(1, 23))

print(f"INFO: Processing {len(ALL_SAMPLES)} samples for WASP analysis")
print(f"INFO: Using 4 VCF groups: {VCF_GROUPS}")

onsuccess:
    print("WASP processing completed successfully!")
    print("Final WASP-processed BAM files: atac_output/wasp/06_blk_filter/*.sorted_final.bam")

rule all:
    input:
        # Chromosome-split VCF files
        expand("/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/{group}/wasp_vcf/chr{chr}.{group}.rename.vcf.gz",
               group=VCF_GROUPS, chr=CHROMOSOMES),
        # HDF5 files for all groups
        expand(OUT("h5", "{group}", "haplotypes.h5"), group=VCF_GROUPS),
        expand(OUT("h5", "{group}", "snp_index.h5"), group=VCF_GROUPS),
        expand(OUT("h5", "{group}", "snp_tab.h5"), group=VCF_GROUPS),
        # WASP filtered BAMs - UPDATED to new step-based directory
        expand(OUT("wasp", "06_blk_filter", "{sampleName}.sorted_final.bam"), sampleName=ALL_SAMPLES),
        expand(OUT("wasp", "06_blk_filter", "{sampleName}.sorted_final.bam.bai"), sampleName=ALL_SAMPLES),
        # QC outputs - UPDATED to new step-based directory  
        expand(OUT("wasp", "07_qc", "ataqv", "{sampleName}.ataqv.json"), sampleName=ALL_SAMPLES),
        expand(OUT("wasp", "07_qc", "bamQC", "{sampleName}_bamQC"), sampleName=ALL_SAMPLES),
        # MACS2 peaks - UPDATED to new step-based directory
        expand(OUT("wasp", "08_peaks", "macs2", "{sampleName}_peaks.narrowPeak"), sampleName=ALL_SAMPLES),
        # Signal tracks - UPDATED to new step-based directory
        expand(OUT("wasp", "09_signals", "indv", "{sampleName}.bw"), sampleName=ALL_SAMPLES)

rule split_vcf_by_chromosome:
    input:
        vcf="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/{group}/{group}.vcf.gz",
        tbi="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/{group}/{group}.vcf.gz.tbi"
    output:
        chr_vcf="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/{group}/wasp_vcf/chr{chr}.{group}.rename.vcf.gz",
        chr_tbi="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/{group}/wasp_vcf/chr{chr}.{group}.rename.vcf.gz.tbi"
    params:
        bcftools_ver=config['bcftools'],
        samtools_version=config["samtoolsVers"],
    log:
        OUT("logs", "split_vcf_{group}_chr{chr}.log")
    shell:
        """
        module load bcftools/{params.bcftools_ver}
        module load samtools/{params.samtools_version}
        mkdir -p $(dirname {output.chr_vcf})
        
        # Extract specific chromosome
        bcftools view -r chr{wildcards.chr} {input.vcf} -O z -o {output.chr_vcf} 2> {log}
        
        # Index chromosome-specific VCF
        tabix -p vcf {output.chr_vcf} 2>> {log}
        """

rule extract_sample_ids:
    input:
        pbs_geno_vcf="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/pbs_geno/pbs_geno.vcf.gz",
        fnf_geno_vcf="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/fnf_geno/fnf_geno.vcf.gz",
        pbs_special_vcf="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/pbs_special/pbs_special.vcf.gz",
        fnf_special_vcf="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/fnf_special/fnf_special.vcf.gz"
    output:
        pbs_geno_samples=OUT("sample_lists", "pbs_geno_sample_ids.txt"),
        fnf_geno_samples=OUT("sample_lists", "fnf_geno_sample_ids.txt"),
        pbs_special_samples=OUT("sample_lists", "pbs_special_sample_ids.txt"),
        fnf_special_samples=OUT("sample_lists", "fnf_special_sample_ids.txt")
    params:
        bcftools_ver=config['bcftools']
    shell:
        """
        module load bcftools/{params.bcftools_ver}
        mkdir -p $(dirname {output.pbs_geno_samples})
        
        # Extract sample names from each group's VCF
        bcftools query -l {input.pbs_geno_vcf} > {output.pbs_geno_samples}
        bcftools query -l {input.fnf_geno_vcf} > {output.fnf_geno_samples}
        bcftools query -l {input.pbs_special_vcf} > {output.pbs_special_samples}
        bcftools query -l {input.fnf_special_vcf} > {output.fnf_special_samples}
        """

rule snp2h5:
    input:
        chr_vcfs=expand("/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/{{group}}/wasp_vcf/chr{chr}.{{group}}.rename.vcf.gz", 
                       chr=CHROMOSOMES)
    output:
        haplotypes=OUT("h5", "{group}", "haplotypes.h5"),
        snp_index=OUT("h5", "{group}", "snp_index.h5"),
        snp_tab=OUT("h5", "{group}", "snp_tab.h5")
    params:
        wasp_path=config['wasp_options']['wasp_dir'],
        chromsize=config['chrsize'],
        wasp_ver=config['waspVer']
    log:
        err=OUT("logs", "{group}_snp2h5.err")
    shell:
        """
        module add wasp/{params.wasp_ver}
        mkdir -p $(dirname {output.haplotypes})

        # Process all chromosome VCF files for this group
        {params.wasp_path}/snp2h5/snp2h5 \
        --chrom {params.chromsize} \
        --format vcf \
        --haplotype {output.haplotypes} \
        --snp_index {output.snp_index} \
        --snp_tab {output.snp_tab} \
        {input.chr_vcfs} 2> {log.err}
        """

rule find_intersecting_snps:
    input:
        bam=lambda wc: BAM_FILES[wc.sampleName],
        bai=lambda wc: BAI_FILES[wc.sampleName],
        # All HDF5 files from all groups
        pbs_geno_snp_tab=OUT("h5", "pbs_geno", "snp_tab.h5"),
        pbs_geno_snp_index=OUT("h5", "pbs_geno", "snp_index.h5"),
        pbs_geno_haplotypes=OUT("h5", "pbs_geno", "haplotypes.h5"),
        fnf_geno_snp_tab=OUT("h5", "fnf_geno", "snp_tab.h5"),
        fnf_geno_snp_index=OUT("h5", "fnf_geno", "snp_index.h5"),
        fnf_geno_haplotypes=OUT("h5", "fnf_geno", "haplotypes.h5"),
        pbs_special_snp_tab=OUT("h5", "pbs_special", "snp_tab.h5"),
        pbs_special_snp_index=OUT("h5", "pbs_special", "snp_index.h5"),
        pbs_special_haplotypes=OUT("h5", "pbs_special", "haplotypes.h5"),
        fnf_special_snp_tab=OUT("h5", "fnf_special", "snp_tab.h5"),
        fnf_special_snp_index=OUT("h5", "fnf_special", "snp_index.h5"),
        fnf_special_haplotypes=OUT("h5", "fnf_special", "haplotypes.h5"),
        # All sample ID files
        pbs_geno_samples=OUT("sample_lists", "pbs_geno_sample_ids.txt"),
        fnf_geno_samples=OUT("sample_lists", "fnf_geno_sample_ids.txt"),
        pbs_special_samples=OUT("sample_lists", "pbs_special_sample_ids.txt"),
        fnf_special_samples=OUT("sample_lists", "fnf_special_sample_ids.txt")
    output:
        # STEP 1: find_intersecting_snps outputs
        remap_fq1=OUT("wasp", "01_find_intersect", "{sampleName}.sorted_final.remap.fq1.gz"),
        remap_fq2=OUT("wasp", "01_find_intersect", "{sampleName}.sorted_final.remap.fq2.gz"),
        to_remap_bam=OUT("wasp", "01_find_intersect", "{sampleName}.sorted_final.to.remap.bam"),
        keep_bam=OUT("wasp", "01_find_intersect", "{sampleName}.sorted_final.keep.bam")
    threads: 4
    params:
        wasp_path=config['wasp_options']['wasp_dir'],
        wasp_ver=config['waspVer'],
        output_dir=OUT("wasp", "01_find_intersect")
    log:
        err=OUT("logs", "01_find_intersecting_snps_{sampleName}.err")
    shell:
        """
        module add wasp/{params.wasp_ver}
        mkdir -p {params.output_dir}

        # Determine which VCF group to use based on sample name characteristics
        if [[ "{wildcards.sampleName}" == *"CTL"* ]] && [[ "{wildcards.sampleName}" != *"replicate"* ]] && [[ "{wildcards.sampleName}" != *"FE"* ]]; then
            # pbs_geno: CTL standard samples
            snp_tab={input.pbs_geno_snp_tab}
            snp_index={input.pbs_geno_snp_index}
            haplotypes={input.pbs_geno_haplotypes}
            sample_file={input.pbs_geno_samples}
            echo "Using pbs_geno for {wildcards.sampleName}"
        elif [[ "{wildcards.sampleName}" == *"FNF"* ]] && [[ "{wildcards.sampleName}" != *"replicate"* ]] && [[ "{wildcards.sampleName}" != *"FE"* ]]; then
            # fnf_geno: FNF standard samples
            snp_tab={input.fnf_geno_snp_tab}
            snp_index={input.fnf_geno_snp_index}
            haplotypes={input.fnf_geno_haplotypes}
            sample_file={input.fnf_geno_samples}
            echo "Using fnf_geno for {wildcards.sampleName}"
        elif [[ "{wildcards.sampleName}" == *"CTL"* ]] && ([[ "{wildcards.sampleName}" == *"replicate"* ]] || [[ "{wildcards.sampleName}" == *"FE"* ]]); then
            # pbs_special: CTL replicates or femur samples
            snp_tab={input.pbs_special_snp_tab}
            snp_index={input.pbs_special_snp_index}
            haplotypes={input.pbs_special_haplotypes}
            sample_file={input.pbs_special_samples}
            echo "Using pbs_special for {wildcards.sampleName}"
        elif [[ "{wildcards.sampleName}" == *"FNF"* ]] && ([[ "{wildcards.sampleName}" == *"replicate"* ]] || [[ "{wildcards.sampleName}" == *"FE"* ]]); then
            # fnf_special: FNF replicates or femur samples
            snp_tab={input.fnf_special_snp_tab}
            snp_index={input.fnf_special_snp_index}
            haplotypes={input.fnf_special_haplotypes}
            sample_file={input.fnf_special_samples}
            echo "Using fnf_special for {wildcards.sampleName}"
        else
            echo "Error: Cannot determine VCF group for sample {wildcards.sampleName}"
            exit 1
        fi

        # Run WASP find_intersecting_snps
        python {params.wasp_path}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir {params.output_dir} \
        --snp_tab $snp_tab \
        --snp_index $snp_index \
        --haplotype $haplotypes \
        --samples $sample_file \
        {input.bam} 2> {log.err}
        """

rule remap_wasp_reads:
    input:
        fq1=rules.find_intersecting_snps.output.remap_fq1,
        fq2=rules.find_intersecting_snps.output.remap_fq2
    output:
        # STEP 2: remap outputs
        sorted_bam=OUT("wasp", "02_remap", "{sampleName}.sorted.bam"),
        bai=OUT("wasp", "02_remap", "{sampleName}.sorted.bam.bai"),
        stats=OUT("wasp", "02_remap", "{sampleName}_stats.txt")
    params:
        bwa_index=config['bwa_index'],
        bwa_version=config['bwaVers'],
        samtools_version=config['samtoolsVers']
    threads: 8
    log:
        err=OUT("logs", "02_remap_wasp_{sampleName}.err")
    shell:
        """
        module load bwa/{params.bwa_version}
        module load samtools/{params.samtools_version}
        mkdir -p $(dirname {output.sorted_bam})
        
        # Remap with BWA (use same parameters as original mapping)
        bwa mem -t {threads} -M {params.bwa_index} {input.fq1} {input.fq2} | \
        samtools view -q 30 -b | \
        samtools sort -o {output.sorted_bam} 2> {log.err}

        # Index and get stats
        samtools index {output.sorted_bam} 2>> {log.err}
        samtools flagstat {output.sorted_bam} > {output.stats} 2>> {log.err}
        """

rule filter_remapped_reads:
    input:
        to_remap_bam=rules.find_intersecting_snps.output.to_remap_bam,
        remapped_bam=rules.remap_wasp_reads.output.sorted_bam
    output:
        keep_bam=OUT("wasp", "03_filter_remapped", "{sampleName}.keep.bam")
    params:
        wasp_path=config['wasp_options']['wasp_dir'],
        wasp_version=config['waspVer']
    log:
        OUT("logs", "03_filter_remapped_{sampleName}.log")
    shell:
        """
        module add wasp/{params.wasp_version}
        mkdir -p $(dirname {output.keep_bam})

        python {params.wasp_path}/mapping/filter_remapped_reads.py \
            {input.to_remap_bam} \
            {input.remapped_bam} \
            {output.keep_bam} 2> {log}
        """

rule merge_keep_bams:
    input:
        remapped_keep=rules.filter_remapped_reads.output.keep_bam,
        original_keep=rules.find_intersecting_snps.output.keep_bam
    output:
        merged_bam=OUT("wasp", "04_merge", "{sampleName}.merged.bam"),
        sorted_bam=OUT("wasp", "04_merge", "{sampleName}.sorted.bam"),
        sorted_bai=OUT("wasp", "04_merge", "{sampleName}.sorted.bam.bai")
    params:
        samtools_version=config['samtoolsVers']
    threads: 4
    log:
        OUT("logs", "04_merge_{sampleName}.log")
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p $(dirname {output.merged_bam})

        # Merge the two keep files
        samtools merge {output.merged_bam} {input.remapped_keep} {input.original_keep} 2> {log}

        # Sort and index
        samtools sort -o {output.sorted_bam} {output.merged_bam} 2>> {log}
        samtools index {output.sorted_bam} 2>> {log}
        """

rule rmdup_wasp:
    input:
        merged_bam=rules.merge_keep_bams.output.sorted_bam
    output:
        rmdup_bam=OUT("wasp", "05_rmdup", "{sampleName}.rmdup.bam"),
        sorted_bam=OUT("wasp", "05_rmdup", "{sampleName}.sorted.bam"),
        sorted_bai=OUT("wasp", "05_rmdup", "{sampleName}.sorted.bam.bai"),
        stats=OUT("wasp", "05_rmdup", "{sampleName}_stats.txt")
    params:
        python_ver=config['python'],
        wasp_path=config['wasp_options']['wasp_dir'],
        wasp_version=config['waspVer'],
        samtools_version=config['samtoolsVers']
    threads: 8
    log:
        err=OUT("logs", "05_rmdup_wasp_{sampleName}.err")
    shell:
        """
        module add wasp/{params.wasp_version}
        module load samtools/{params.samtools_version}
        module load python/{params.python_ver}
        mkdir -p $(dirname {output.rmdup_bam})

        # Use WASP's unbiased duplicate removal
        python {params.wasp_path}/mapping/rmdup_pe.py {input.merged_bam} {output.rmdup_bam} 2> {log.err}

        # Sort and index
        samtools sort -o {output.sorted_bam} {output.rmdup_bam} 2>> {log.err}
        samtools index {output.sorted_bam} 2>> {log.err}
        samtools flagstat {output.sorted_bam} > {output.stats} 2>> {log.err}
        """

rule remove_blacklist_wasp:
    input:
        bam=rules.rmdup_wasp.output.sorted_bam
    output:
        final_bam=OUT("wasp", "06_blk_filter", "{sampleName}.sorted_final.bam"),
        final_bai=OUT("wasp", "06_blk_filter", "{sampleName}.sorted_final.bam.bai"),
        stats=OUT("wasp", "06_blk_filter", "{sampleName}_stats.txt")
    threads: 4
    params:
        bedtools_ver=config['bedtools'],
        samtools_version=config['samtoolsVers'],
        blacklist_bed=config['blacklist']
    log:
        err=OUT("logs", "06_blacklist_filter_wasp_{sampleName}.err")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        module load samtools/{params.samtools_version}
        mkdir -p $(dirname {output.final_bam})

        # Remove blacklist regions
        bedtools intersect -v -abam {input.bam} -b {params.blacklist_bed} | \
          samtools sort -o {output.final_bam} 2> {log.err}

        # Index and get stats
        samtools index {output.final_bam} 2>> {log.err}
        samtools flagstat {output.final_bam} > {output.stats} 2>> {log.err}
        """

rule run_ataqv_wasp:
    input:
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        txt=OUT("wasp", "07_qc", "ataqv", "{sampleName}.ataqv.txt"),
        json=OUT("wasp", "07_qc", "ataqv", "{sampleName}.ataqv.json")
    params:
        tss_file=config["tss_annotation"],
        ataqv_bin=config["ataqvVers"],
        tss_extension=config["tssextension"]
    threads: 6
    log:
        err=OUT("logs", "07_ataqv_wasp_{sampleName}.err")
    shell:
        """
        mkdir -p $(dirname {output.txt})

        {params.ataqv_bin} --tss-file {params.tss_file} \
            --tss-extension {params.tss_extension} \
            --metrics-file {output.json} \
            --ignore-read-groups human \
            {input.bam} > {output.txt} 2> {log.err}
        """

rule bam_qc_wasp:
    input:
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        qc=OUT("wasp", "07_qc", "bamQC", "{sampleName}_bamQC")
    params:
        script=config["atac_qc_Rscript"],
        r_version=config["r"],
        output_dir=OUT("wasp", "07_qc", "bamQC")
    threads: 6
    log:
        err=OUT("logs", "07_bam_qc_wasp_{sampleName}.err")
    shell:
        """
        module load r/{params.r_version}
        mkdir -p {params.output_dir}

        Rscript {params.script} {input.bam} "{params.output_dir}" ".sorted_final.bam" 2> {log.err}
        """

rule macs2_peaks:
    input:
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        peak=OUT("wasp", "08_peaks", "macs2", "{sampleName}_peaks.narrowPeak")
    log:
        err=OUT("logs", "08_macs2_{sampleName}.err")
    threads: 2
    params:
        python_ver=config['python'],
        macs2_ver=config['macs2'],
        output_dir=OUT("wasp", "08_peaks", "macs2")
    shell:
        """
        module load python/{params.python_ver}
        module load macs/{params.macs2_ver}
        mkdir -p {params.output_dir}

        macs2 callpeak -t {input.bam} \
            -f BAMPE \
            -g hs \
            --nomodel \
            --shift -75 \
            --extsize 150 \
            --keep-dup all \
            -p 0.01 \
            --call-summits \
            -n {wildcards.sampleName} \
            --outdir {params.output_dir} 2> {log.err}
        """

rule signal_tracks:
    input:
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        signal=OUT("wasp", "09_signals", "indv", "{sampleName}.bw")
    log:
        err=OUT("logs", "09_signal_{sampleName}.err")
    params:
        deeptools_ver=config['deeptoolsVers'],
        bin_size=config['bin_size'],
        effective_genome_size=config['effective_genome_size'],
        normalize_option=config['normalize_option']
    threads: 4
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p $(dirname {output.signal})

        bamCoverage \
            --bam {input.bam} \
            -o {output.signal} \
            --binSize {params.bin_size} \
            --normalizeUsing {params.normalize_option} \
            --effectiveGenomeSize {params.effective_genome_size} \
            --extendReads \
            --minMappingQuality 30 \
            --ignoreForNormalization chrX chrY chrM \
            --numberOfProcessors {threads} 2> {log.err}
        """
