#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

#------------------------------------------------------
# Load sample sheet (same logic as your WASP pipeline)
bam_samplesheet_path = os.path.join(config["paths"]["outdir"], config["bam_samplesheet"])
bam_samples = pd.read_csv(bam_samplesheet_path, sep='\t')
samples_df = bam_samples.copy()
s_config = config['sample_logic']
donor_col = s_config['femur_donor_col']

#------------------------------------------------------
# Femur or replicate samples
if str(s_config.get('Sample_handing', 'false')).lower() == 'true':
    rep_col = s_config['replicate_column']
    is_replicate = samples_df[rep_col].str.contains('replicate', case=False, na=False)
    
    tissue_col = s_config['femur_column']
    is_femur = samples_df[tissue_col].str.contains('Femur', na=False) | \
               samples_df[donor_col].str.endswith('FE', na=False)
    
    # Full name (used in BAMs)
    full_df = samples_df.copy()
    full_df['replicate_suffix'] = np.where(is_replicate, 'replicate', '1')
    full_df[donor_col] = np.where(
        is_femur & ~full_df[donor_col].str.endswith('FE'),
        full_df[donor_col] + 'FE',
        full_df[donor_col]
    )
    full_name_cols = s_config['base_name_columns'] + ['replicate_suffix']
    samples_df['full_sample_name'] = full_df[full_name_cols].astype(str).agg('_'.join, axis=1)
    
    is_special = is_replicate | is_femur
else:
    std_name = samples_df[s_config['base_name_columns']].astype(str).agg('_'.join, axis=1) + '_1'
    samples_df['full_sample_name'] = std_name
    is_special = pd.Series(False, index=samples_df.index)

group_col = s_config['condition_column']
ALL_SAMPLES = samples_df['full_sample_name'].tolist()

# Getting the bam
BAM_FILES = {}
BAI_FILES = {}
for _, row in samples_df.iterrows():
    if pd.notna(row.get('Bam_directory')) and pd.notna(row.get('Bam_file')):
        sample_name = row['full_sample_name']
        BAM_FILES[sample_name] = os.path.join(row['Bam_directory'], row['Bam_file'])
        BAI_FILES[sample_name] = BAM_FILES[sample_name] + ".bai"


#for _, row in samples_df.iterrows():
#    sample_name = row['full_sample_name']
#    bam_path = OUT("filtered", "blk_filter", f"{sample_name}.sorted_final.bam")
#    BAM_FILES[sample_name] = bam_path
#    BAI_FILES[sample_name] = BAM_FILES[sample_name] + ".bai"

MAIN_SAMPLES = samples_df[~is_special]['full_sample_name'].tolist()
MAIN_PBS = samples_df[(samples_df[group_col] == s_config['control_group_value']) & ~is_special]['full_sample_name'].tolist()
MAIN_FNF = samples_df[(samples_df[group_col] == s_config['treatment_group_value']) & ~is_special]['full_sample_name'].tolist()

print(f"INFO: {len(ALL_SAMPLES)} total samples")
print(f"INFO: {len(MAIN_SAMPLES)} main samples (21 donors)")
print(f"INFO: PBS: {len(MAIN_PBS)}, FNF: {len(MAIN_FNF)}")

#------------------------------------------------------
rule all:
    input:
        # 1. Individual peak calling for ALL samples
        expand(OUT("peaks", "individual", "{sampleName}_peaks.narrowPeak"),
               sampleName=ALL_SAMPLES),

        # 2. Individual peak counts for ALL samples
        expand(OUT("peaks", "counts", "individual", "{sampleName}_peaks_counts.txt"),
               sampleName=ALL_SAMPLES),

        # 3. Merged peaks for PBS and FNF separately (no femur or replicate)
        OUT("peaks", "merged", "PBS_merged.bed"),
        OUT("peaks", "merged", "FNF_merged.bed"),

        # 4. Peak counts on merged peaks for PBS and FNF separately
        OUT("peaks", "counts", "PBS_merged_counts.txt"),
        OUT("peaks", "counts", "FNF_merged_counts.txt"),

        # 5. All samples big peak count (merged from condition)
        OUT("peaks", "counts", "all_peaks_combined_counts.txt"),

        # 6. FRiP scores for ALL samples
        expand(OUT("peaks", "frip", "{sampleName}_frip.txt"),
               sampleName=ALL_SAMPLES),
        expand(OUT("peaks", "frip", "{sampleName}_peaks_counts.txt"), sampleName=ALL_SAMPLES),
        # 7. Individual signal tracks for ALL samples
        expand(OUT("signals", "individual", "{sampleName}.bw"),
               sampleName=ALL_SAMPLES),

        # 8. Merged signal tracks for PBS and FNF (no femur or replicate)
        OUT("signals", "merged", "PBS_merged.bw"),
        OUT("signals", "merged", "FNF_merged.bw")
#------------------------------------------------------
rule call_indiv_peaks:
    input:
        bam=lambda wc: BAM_FILES[wc.sampleName]
    output:
        peak=OUT("peaks", "individual", "{sampleName}_peaks.narrowPeak"),
        summit=OUT("peaks", "individual", "{sampleName}_summits.bed")
    log:
        err=OUT("logs", "peaks_{sampleName}.err")
    threads: 2
    params:
        python_ver=config['python'],
        macs2_ver=config['macs2'],
        output_dir=OUT("peaks", "individual")
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

rule count_indiv_peaks:
    input:
        peaks=rules.call_indiv_peaks.output.peak,
        bam=lambda wc: BAM_FILES[wc.sampleName]
    output:
        saf=OUT("peaks", "counts", "individual", "{sampleName}_peaks.saf"),
        counts=OUT("peaks", "counts", "individual", "{sampleName}_peaks_counts.txt")
    params:
        subread_ver=config['subread']
    threads: 4
    log:
        err=OUT("logs", "count_individual_{sampleName}.err")
    shell:
        """
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.saf})

        # Convert peaks to SAF format
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {input.peaks} > {output.saf}

        # Count reads in individual peaks
        featureCounts -p -F SAF -T {threads} --primary -C \
            -a {output.saf} -o {output.counts} {input.bam} 2> {log.err}
        """

rule merge_peaks_PBS:
    input:
        peaks=expand(OUT("peaks", "individual", "{sample}_peaks.narrowPeak"), sample=MAIN_PBS)
    output:
        bed=OUT("peaks", "merged", "PBS_merged.bed"),
        saf=OUT("peaks", "merged", "PBS_merged.saf")
    params:
        bedtools_ver=config['bedtools']
    log:
        err=OUT("logs", "merge_PBS.err")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        mkdir -p $(dirname {output.bed})

        # Merge all PBS peaks (union), filter for size
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge | \
        awk '($3-$2) >= 50 && ($3-$2) <= 3000' > {output.bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
             {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}
        """

rule merge_peaks_FNF:
    input:
        peaks=expand(OUT("peaks", "individual", "{sample}_peaks.narrowPeak"), sample=MAIN_FNF)
    output:
        bed=OUT("peaks", "merged", "FNF_merged.bed"),
        saf=OUT("peaks", "merged", "FNF_merged.saf")
    params:
        bedtools_ver=config['bedtools']
    log:
        err=OUT("logs", "merge_FNF.err")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        mkdir -p $(dirname {output.bed})

        # Merge all FNF peaks (union), filter for size
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge | \
        awk '($3-$2) >= 50 && ($3-$2) <= 3000' > {output.bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}
        """

rule count_PBS_merged_peaks:
    input:
        saf=rules.merge_peaks_PBS.output.saf,
        bams=[BAM_FILES[sample] for sample in MAIN_PBS]
    output:
        counts=OUT("peaks", "counts", "PBS_merged_counts.txt")
    params:
        subread_ver=config['subread']
    threads: 8
    log:
        err=OUT("logs", "count_PBS_merged.err")
    shell:
        """
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.counts})

        featureCounts -p -F SAF -T {threads} --primary -C \
            -a {input.saf} -o {output.counts} {input.bams} 2> {log.err}
        """

rule count_FNF_merged_peaks:
    input:
        saf=rules.merge_peaks_FNF.output.saf,
        bams=[BAM_FILES[sample] for sample in MAIN_FNF]
    output:
        counts=OUT("peaks", "counts", "FNF_merged_counts.txt")
    params:
        subread_ver=config['subread']
    threads: 8
    log:
        err=OUT("logs", "count_FNF_merged.err")
    shell:
        """
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.counts})

        featureCounts -p -F SAF -T {threads} --primary -C \
            -a {input.saf} -o {output.counts} {input.bams} 2> {log.err}
        """

rule peak_count_all:
    input:
        pbs_bed=rules.merge_peaks_PBS.output.bed,
        fnf_bed=rules.merge_peaks_FNF.output.bed,
        bams=[BAM_FILES[sample] for sample in ALL_SAMPLES]
    output:
        big_bed=OUT("peaks", "merged", "all_peaks_combined.bed"),
        big_saf=OUT("peaks", "merged", "all_peaks_combined.saf"),
        big_counts=OUT("peaks", "counts", "all_peaks_combined_counts.txt")
    params:
        bedtools_ver=config['bedtools'],
        subread_ver=config['subread']
    threads: 8
    log:
        err=OUT("logs", "big_peaks_count.err")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.big_bed})

        # Combine PBS and FNF peaks and merge overlaps
        cat {input.pbs_bed} {input.fnf_bed} | \
        sort -k1,1 -k2,2n | \
        bedtools merge > {output.big_bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.big_bed} > {output.big_saf}

        # Count ALL samples on big peak set
        featureCounts -p -F SAF -T {threads} --primary -C \
            -a {output.big_saf} -o {output.big_counts} {input.bams} 2>> {log.err}
        """

rule calculate_frip_all:
    input:
        saf=rules.peak_count_all.output.big_saf,
        bam=lambda wc: BAM_FILES[wc.sampleName]
    output:
        frip=OUT("peaks", "frip", "{sampleName}_frip.txt"),
        counts=OUT("peaks", "frip", "{sampleName}_peaks_counts.txt")  # Keep the counts too
    params:
        subread_ver=config['subread']
    threads: 4
    log:
        err=OUT("logs", "frip_{sampleName}.err")
    shell:
        """
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.frip})

        # Run featureCounts - it will create .summary file automatically
        featureCounts -p -F SAF -T {threads} --primary -C --fracOverlap 0.2 \
            -a {input.saf} -o {output.counts} {input.bam} 2> {log.err}

        # Extract FRiP from summary file
        summary_file="{output.counts}.summary"
        assigned=$(awk 'NR==2 {{print $2}}' $summary_file)
        no_features=$(awk '/Unassigned_NoFeatures/ {{print $2}}' $summary_file)
        overlapping=$(awk '/Unassigned_Overlapping_Length/ {{print $2}}' $summary_file)
        
        total_relevant=$((assigned + no_features + overlapping))
        frip=$(echo "scale=4; $assigned / $total_relevant" | bc -l)

        echo -e "sample\tassigned\tno_features\toverlapping\ttotal_relevant\tfrip_score" > {output.frip}
        echo -e "{wildcards.sampleName}\t$assigned\t$no_features\t$overlapping\t$total_relevant\t$frip" >> {output.frip}
        """

rule create_individual_signals:
    input:
        bam=lambda wc: BAM_FILES[wc.sampleName]
    output:
        signal=OUT("signals", "individual", "{sampleName}.bw")
    log:
        err=OUT("logs", "signal_{sampleName}.err")
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

rule merge_bams_PBS:
    input:
        bams=[BAM_FILES[sample] for sample in MAIN_PBS]
    output:
        merged_bam=OUT("signals", "merged", "PBS_merged.bam"),
        merged_bai=OUT("signals", "merged", "PBS_merged.bam.bai")
    params:
        samtools_ver=config['samtoolsVers']
    threads: 8
    log:
        err=OUT("logs", "merge_bam_PBS.err")
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p $(dirname {output.merged_bam})

        samtools merge -@ {threads} {output.merged_bam} {input.bams} 2> {log.err}
        samtools index {output.merged_bam} 2>> {log.err}
        """

rule merge_bams_FNF:
    input:
        bams=[BAM_FILES[sample] for sample in MAIN_FNF]
    output:
        merged_bam=OUT("signals", "merged", "FNF_merged.bam"),
        merged_bai=OUT("signals", "merged", "FNF_merged.bam.bai")
    params:
        samtools_ver=config['samtoolsVers']
    threads: 8
    log:
        err=OUT("logs", "merge_bam_FNF.err")
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p $(dirname {output.merged_bam})

        samtools merge -@ {threads} {output.merged_bam} {input.bams} 2> {log.err}
        samtools index {output.merged_bam} 2>> {log.err}
        """

rule merged_signals:
    input:
        bam=OUT("signals", "merged", "{condition}_merged.bam")
    output:
        signal=OUT("signals", "merged", "{condition}_merged.bw")
    log:
        err=OUT("logs", "merged_signal_{condition}.err")
    params:
        deeptools_ver=config['deeptoolsVers'],
        bin_size=config['bin_size'],
        effective_genome_size=config['effective_genome_size'],
        normalize_option=config['normalize_option']
    threads: 8
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
