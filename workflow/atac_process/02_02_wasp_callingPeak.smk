#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

# Include the WASP workflow to reuse all sample processing logic
include: "02_wasp.smk"

# All variables from 02_wasp.smk are now available:
# - ALL_SAMPLES, BAM_FILES, BAI_FILES, samples, is_replicate, is_femur, etc.

# We need to define the groupings for peak merging that aren't in 02_wasp.smk
is_special = is_replicate | is_femur if str(s_config.get('Sample_handing', 'false')).lower() == 'true' else pd.Series(False, index=samples.index)
MAIN_SAMPLES = samples[~is_special]['full_sample_name'].tolist()
MAIN_PBS = samples[(samples[group_col] == s_config['control_group_value']) & ~is_special]['full_sample_name'].tolist()
MAIN_FNF = samples[(samples[group_col] == s_config['treatment_group_value']) & ~is_special]['full_sample_name'].tolist()

print(f"INFO: WASP Peak Calling - Processing {len(ALL_SAMPLES)} total samples")
print(f"INFO: {len(MAIN_SAMPLES)} main samples (standard donors)")
print(f"INFO: PBS: {len(MAIN_PBS)}, FNF: {len(MAIN_FNF)}")
print(f"INFO: Special samples (femur/replicate): {len(ALL_SAMPLES) - len(MAIN_SAMPLES)}")

#------------------------------------------------------
rule wasp_peak_all:
    input:
        # 1. Individual peak calling for ALL samples (WASP-processed)
        expand(OUT("wasp", "peaks", "individual", "{sampleName}_peaks.narrowPeak"),
               sampleName=ALL_SAMPLES),

        # 2. Individual peak counts for ALL samples
        expand(OUT("wasp", "peaks", "counts", "individual", "{sampleName}_peaks_counts.txt"),
               sampleName=ALL_SAMPLES),

        # 3. Merged peaks for PBS and FNF separately (main samples only - no femur or replicate)
        OUT("wasp", "peaks", "merged", "PBS_merged.bed"),
        OUT("wasp", "peaks", "merged", "FNF_merged.bed"),

        # 4. Peak counts on merged peaks for PBS and FNF separately (main samples only)
        OUT("wasp", "peaks", "counts", "PBS_merged_counts.txt"),
        OUT("wasp", "peaks", "counts", "FNF_merged_counts.txt"),

        # 5. All samples big peak count (merged from all conditions)
        OUT("wasp", "peaks", "counts", "all_peaks_combined_counts.txt"),

        # 6. FRiP scores for ALL samples (including femur and replicates)
        expand(OUT("wasp", "peaks", "frip", "{sampleName}_frip.txt"),
               sampleName=ALL_SAMPLES),
        expand(OUT("wasp", "peaks", "frip", "{sampleName}_peaks_counts.txt"), sampleName=ALL_SAMPLES),
        
        # 7. Individual signal tracks for ALL samples (including femur and replicates)
        expand(OUT("wasp", "signals", "individual", "{sampleName}.bw"),
               sampleName=ALL_SAMPLES),

        # 8. Merged signal tracks for PBS and FNF (main samples only - no femur or replicate)
        OUT("wasp", "signals", "merged", "PBS_merged.bw"),
        OUT("wasp", "signals", "merged", "FNF_merged.bw")

#------------------------------------------------------
rule call_wasp_indiv_peaks:
    input:
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        peak=OUT("wasp", "peaks", "individual", "{sampleName}_peaks.narrowPeak"),
        summit=OUT("wasp", "peaks", "individual", "{sampleName}_summits.bed")
    log:
        err=OUT("logs", "wasp_peaks_{sampleName}.err")
    threads: 2
    params:
        python_ver=config['python'],
        macs2_ver=config['macs2'],
        output_dir=OUT("wasp", "peaks", "individual")
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

rule count_wasp_indiv_peaks:
    input:
        peaks=rules.call_wasp_indiv_peaks.output.peak,
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        saf=OUT("wasp", "peaks", "counts", "individual", "{sampleName}_peaks.saf"),
        counts=OUT("wasp", "peaks", "counts", "individual", "{sampleName}_peaks_counts.txt")
    params:
        subread_ver=config['subread']
    threads: 4
    log:
        err=OUT("logs", "wasp_count_individual_{sampleName}.err")
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

rule merge_wasp_peaks_PBS:
    input:
        peaks=expand(OUT("wasp", "peaks", "individual", "{sample}_peaks.narrowPeak"), sample=MAIN_PBS)
    output:
        bed=OUT("wasp", "peaks", "merged", "PBS_merged.bed"),
        saf=OUT("wasp", "peaks", "merged", "PBS_merged.saf")
    params:
        bedtools_ver=config['bedtools']
    log:
        err=OUT("logs", "wasp_merge_PBS.err")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        mkdir -p $(dirname {output.bed})

        # Use bedtools merge 
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -d 100 | \
        awk '($3-$2) >= 40 && ($3-$2) <= 3000' > {output.bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
             {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}
        """

rule merge_wasp_peaks_FNF:
    input:
        peaks=expand(OUT("wasp", "peaks", "individual", "{sample}_peaks.narrowPeak"), sample=MAIN_FNF)
    output:
        bed=OUT("wasp", "peaks", "merged", "FNF_merged.bed"),
        saf=OUT("wasp", "peaks", "merged", "FNF_merged.saf")
    params:
        bedtools_ver=config['bedtools']
    log:
        err=OUT("logs", "wasp_merge_FNF.err")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        mkdir -p $(dirname {output.bed})

        # Use bedtools merge 
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -d 100 | \
        awk '($3-$2) >= 40 && ($3-$2) <= 3000' > {output.bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}
        """

rule count_wasp_PBS_merged_peaks:
    input:
        saf=rules.merge_wasp_peaks_PBS.output.saf,
        bams=expand(OUT("wasp", "06_blk_filter", "{sample}.sorted_final.bam"), sample=MAIN_PBS)
    output:
        counts=OUT("wasp", "peaks", "counts", "PBS_merged_counts.txt")
    params:
        subread_ver=config['subread']
    threads: 8
    log:
        err=OUT("logs", "wasp_count_PBS_merged.err")
    shell:
        """
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.counts})

        featureCounts -p -F SAF -T {threads} --primary -C \
            -a {input.saf} -o {output.counts} {input.bams} 2> {log.err}
        """

rule count_wasp_FNF_merged_peaks:
    input:
        saf=rules.merge_wasp_peaks_FNF.output.saf,
        bams=expand(OUT("wasp", "06_blk_filter", "{sample}.sorted_final.bam"), sample=MAIN_FNF)
    output:
        counts=OUT("wasp", "peaks", "counts", "FNF_merged_counts.txt")
    params:
        subread_ver=config['subread']
    threads: 8
    log:
        err=OUT("logs", "wasp_count_FNF_merged.err")
    shell:
        """
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.counts})

        featureCounts -p -F SAF -T {threads} --primary -C \
            -a {input.saf} -o {output.counts} {input.bams} 2> {log.err}
        """

rule wasp_peak_count_all:
    input:
        peaks=expand(OUT("wasp", "peaks", "individual", "{sample}_peaks.narrowPeak"), sample=MAIN_SAMPLES),
        bams=expand(OUT("wasp", "06_blk_filter", "{sample}.sorted_final.bam"), sample=MAIN_SAMPLES)
    output:
        big_bed=OUT("wasp", "peaks", "merged", "all_peaks_combined.bed"),
        big_saf=OUT("wasp", "peaks", "merged", "all_peaks_combined.saf"),
        big_counts=OUT("wasp", "peaks", "counts", "all_peaks_combined_counts.txt")
    params:
        bedtools_ver=config['bedtools'],
        subread_ver=config['subread']
    threads: 8
    log:
        err=OUT("logs", "wasp_big_peaks_count.err")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        module load subread/{params.subread_ver}
        mkdir -p $(dirname {output.big_bed})

        # Use bedtools merge 
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -d 100 | \
        awk '($3-$2) >= 40 && ($3-$2) <= 3000' > {output.big_bed} 2> {log.err}
        
        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.big_bed} > {output.big_saf}

        # Count ALL samples on peak set
        featureCounts -p -F SAF -T {threads} --primary -C \
            -a {output.big_saf} -o {output.big_counts} {input.bams} 2>> {log.err}
        """

rule calculate_wasp_frip_all:
    input:
        saf=rules.wasp_peak_count_all.output.big_saf,
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        frip=OUT("wasp", "peaks", "frip", "{sampleName}_frip.txt"),
        counts=OUT("wasp", "peaks", "frip", "{sampleName}_peaks_counts.txt")
    params:
        subread_ver=config['subread']
    threads: 4
    log:
        err=OUT("logs", "wasp_frip_{sampleName}.err")
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

rule create_wasp_individual_signals:
    input:
        bam=rules.remove_blacklist_wasp.output.final_bam
    output:
        signal=OUT("wasp", "signals", "individual", "{sampleName}.bw")
    log:
        err=OUT("logs", "wasp_signal_{sampleName}.err")
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

rule merge_wasp_bams_PBS:
    input:
        bams=expand(OUT("wasp", "06_blk_filter", "{sample}.sorted_final.bam"), sample=MAIN_PBS)
    output:
        merged_bam=OUT("wasp", "signals", "merged", "PBS_merged.bam"),
        merged_bai=OUT("wasp", "signals", "merged", "PBS_merged.bam.bai")
    params:
        samtools_ver=config['samtoolsVers']
    threads: 8
    log:
        err=OUT("logs", "wasp_merge_bam_PBS.err")
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p $(dirname {output.merged_bam})

        samtools merge -@ {threads} {output.merged_bam} {input.bams} 2> {log.err}
        samtools index {output.merged_bam} 2>> {log.err}
        """

rule merge_wasp_bams_FNF:
    input:
        bams=expand(OUT("wasp", "06_blk_filter", "{sample}.sorted_final.bam"), sample=MAIN_FNF)
    output:
        merged_bam=OUT("wasp", "signals", "merged", "FNF_merged.bam"),
        merged_bai=OUT("wasp", "signals", "merged", "FNF_merged.bam.bai")
    params:
        samtools_ver=config['samtoolsVers']
    threads: 8
    log:
        err=OUT("logs", "wasp_merge_bam_FNF.err")
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p $(dirname {output.merged_bam})

        samtools merge -@ {threads} {output.merged_bam} {input.bams} 2> {log.err}
        samtools index {output.merged_bam} 2>> {log.err}
        """

rule wasp_merged_signals:
    input:
        bam=OUT("wasp", "signals", "merged", "{condition}_merged.bam")
    output:
        signal=OUT("wasp", "signals", "merged", "{condition}_merged.bw")
    log:
        err=OUT("logs", "wasp_merged_signal_{condition}.err")
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
