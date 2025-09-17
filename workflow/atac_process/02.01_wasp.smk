#!/usr/bin/env python3

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakefiles.utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
samples = samples.astype(str)

# Exclusion criteria
if config['exclude_options']['exclude_replicates']:
    samples = samples[~samples['Protocol_notes'].str.contains('replicate', case=False, na=False)]
if config['exclude_options']['exclude_synovium']:
    samples = samples[~samples['Tissue'].str.contains('synovium', case=False, na=False)]
if config['exclude_options']['exclude_femur']:
    samples = samples[~samples['Tissue'].str.contains('Femur', case=False, na=False)]

## Create merged names
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ['Proj','Donor','Condition','Tissue','Protocol_notes']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

# Separate samples by condition
#ctl_samples = samples[samples['Condition'] == 'CTL']
#fnf_samples = samples[samples['Condition'] == 'FNF']

# Group samples by condition
condition_samples = samples.groupby('Condition')['mn'].apply(list).to_dict()

# Create BAM and BAI file mappings
bam_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/output/filtered/{x.name}_filtered.bam").to_dict()
bai_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/output/filtered/{x.name}_filtered.bam.bai").to_dict()

rule all:
    input:
        expand("output/wasp/filter_remapped_reads/{sampleName}_filtered.{ext}",
               sampleName=bam_files.keys(),
               ext=['remap.fq1.gz', 'remap.fq2.gz', 'to.remap.bam', 'keep.bam']),
        expand("output/wasp/filtered_bam/{sampleName}.keep.bam",
               sampleName=bam_files.keys()),
        expand("output/wasp/rmdup/{sampleName}{ext}",
               sampleName=bam_files.keys(),
               ext=['.rmdup.bam', '.sorted_wasp.bam', '.sorted_wasp.bam.bai', '_stats.txt']),
        expand("output/wasp/blk_filter/{sampleName}{ext}",
               sampleName=bam_files.keys(),
               ext=[".blkfilter.bam", ".sorted_final.bam", ".sorted_final.bam.bai", "_stats.txt"]),
        expand("output/wasp/ataqv/{sampleName}.ataqv.json", sampleName=bam_files.keys()),
        expand("output/wasp/ataqv/{sampleName}.ataqv.txt", sampleName=bam_files.keys()),
        [expand("output/wasp/bamQC/{sampleName}_bamqQC", sampleName=key) for key in bam_files.keys()],
        expand("output/peaks/merged_cond/{cond}_macs3_merged{ext}",
               ext=['.bed', '.saf', '_counts.txt'], cond=['CTL', 'FNF']),
        expand("output/peaks/merged/allsamples_{tool}_merged.bed", tool=['macs3', 'macs2', 'hmmratac', 'rocco']),
        expand("output/peaks/merged/allsamples_{tool}_merged.saf", tool=['macs3', 'macs2', 'hmmratac', 'rocco']),
        expand("output/peaks/merged/allsamples_{tool}_merged_counts.txt", tool=['macs3', 'macs2', 'hmmratac', 'rocco']),
        expand("output/peaks/{tool}/indv/{sampleName}_count.txt",
               tool=['macs3', 'macs2', 'hmmratac', 'rocco'], sampleName=bam_files.keys()),
        expand("output/signals/merged_signal/{cond}_merged.bw", cond=['CTL', 'FNF']),
        expand("output/signals/indv/{sampleName}.bw", sampleName=bam_files.keys())


rule extract_sample_ids:
    input:
        pbs_matched="output/geno/pbs_geno/pbs_matched.txt",
        fnf_matched="output/geno/fnf_geno/fnf_matched.txt"
    output:
        sample_ids_pbs="output/wasp/pbs_sample_ids.txt",
        sample_ids_fnf="output/wasp/fnf_sample_ids.txt"
    log:
        pbs_log="output/logs/extract_sample_ids_pbs.log",
        fnf_log="output/logs/extract_sample_ids_fnf.log"
    shell:
        """

        # Extract the first column for PBS sample IDs
        cut -d' ' -f1 {input.pbs_matched} > {output.sample_ids_pbs} 2> {log.pbs_log}

        # Extract the first column for FNF sample IDs
        cut -d' ' -f1 {input.fnf_matched} > {output.sample_ids_fnf} 2> {log.fnf_log}
        """

rule snp2h5:
    input:
        sample_ids_pbs=rules.extract_sample_ids.output.sample_ids_pbs,
        sample_ids_fnf=rules.extract_sample_ids.output.sample_ids_fnf
    output:
        pbs_haplotypes="output/wasp/pbs/haplotypes.h5",
        pbs_snp_index="output/wasp/pbs/snp_index.h5",
        pbs_snp_tab="output/wasp/pbs/snp_tab.h5",
        fnf_haplotypes="output/wasp/fnf/haplotypes.h5",
        fnf_snp_index="output/wasp/fnf/snp_index.h5",
        fnf_snp_tab="output/wasp/fnf/snp_tab.h5"
    params:
        wasp_path=config['wasp_dir'],
        chromsize=config['chrom_size'],
        wasp_ver=config['waspVer']
    log:
        pbs_err="output/logs/pbs_snp2h5.err",
        fnf_err="output/logs/fnf_snp2h5.err"
    shell:
        """
        module add wasp/{params.wasp_ver}
        mkdir -p output/wasp/pbs output/wasp/fnf

        # Process PBS VCF
        {params.wasp_path}/snp2h5/snp2h5 \
        --chrom {params.chromsize} \
        --format vcf \
        --haplotype {output.pbs_haplotypes} \
        --snp_index {output.pbs_snp_index} \
        --snp_tab {output.pbs_snp_tab} \
        output/geno/pbs_geno/wasp_vcf/chr*.pbs.rename.vcf.gz 2> {log.pbs_err}

        # Process FNF VCF
        {params.wasp_path}/snp2h5/snp2h5 \
        --chrom {params.chromsize} \
        --format vcf \
        --haplotype {output.fnf_haplotypes} \
        --snp_index {output.fnf_snp_index} \
        --snp_tab {output.fnf_snp_tab} \
        output/geno/fnf_geno/wasp_vcf/chr*.fnf.rename.vcf.gz 2> {log.fnf_err}
        """

rule find_intersecting_snps:
    input:
        bam=lambda wildcards: bam_files[wildcards.sampleName],
        pbs_snp_tab=rules.snp2h5.output.pbs_snp_tab,
        pbs_snp_index=rules.snp2h5.output.pbs_snp_index,
        pbs_haplotypes=rules.snp2h5.output.pbs_haplotypes,
        fnf_snp_tab=rules.snp2h5.output.fnf_snp_tab,
        fnf_snp_index=rules.snp2h5.output.fnf_snp_index,
        fnf_haplotypes=rules.snp2h5.output.fnf_haplotypes,
        sample_ids_pbs=rules.extract_sample_ids.output.sample_ids_pbs,
        sample_ids_fnf=rules.extract_sample_ids.output.sample_ids_fnf
    output:
        pbs_fq1="output/wasp/filter_remapped_reads/{sampleName}_filtered.remap.fq1.gz",
        pbs_fq2="output/wasp/filter_remapped_reads/{sampleName}_filtered.remap.fq2.gz",
        to_remap_bam="output/wasp/filter_remapped_reads/{sampleName}_filtered.to.remap.bam",
        keep_bam="output/wasp/filter_remapped_reads/{sampleName}_filtered.keep.bam"
    threads: 4
    params:
        wasp_path=config['wasp_dir'],
        wasp_ver=config['waspVer']
    log:
        err="output/logs/find_intersecting_snps_pbs_{sampleName}.err"
    shell:
        """
        module add wasp/{params.wasp_ver}

        mkdir -p output/wasp/filter_remapped_reads

        if [[ "{wildcards.sampleName}" == *"CTL"* ]]; then
            snp_tab={input.pbs_snp_tab}
            snp_index={input.pbs_snp_index}
            haplotypes={input.pbs_haplotypes}
            sample_id={input.sample_ids_pbs}
        elif [[ "{wildcards.sampleName}" == *"FNF"* ]]; then
            snp_tab={input.fnf_snp_tab}
            snp_index={input.fnf_snp_index}
            haplotypes={input.fnf_haplotypes}
            sample_id={input.sample_ids_fnf}
        else
            echo "Error: Sample type not recognized in {wildcards.sampleName}"
            exit 1
        fi

        # Add error checking here
        if [[ -z "$snp_tab" ]]; then
            echo "Error: snp_tab not set"
            exit 1
        fi
        if [[ -z "$snp_index" ]]; then
            echo "Error: snp_index not set"
            exit 1
        fi
        if [[ -z "$haplotypes" ]]; then
            echo "Error: haplotypes not set"
            exit 1
        fi
        if [[ -z "$sample_id" ]]; then
            echo "Error: sample_id not set"
            exit 1
        fi 

        # Then run WASP command
        python {params.wasp_path}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir output/wasp/filter_remapped_reads \
        --snp_tab $snp_tab \
        --snp_index $snp_index \
        --haplotype $haplotypes \
        --samples $sample_id \
        {input.bam} 2> {log.err}

        """


rule remap_wasp_reads:
    input:
        fq1=rules.find_intersecting_snps.output.pbs_fq1,
        fq2=rules.find_intersecting_snps.output.pbs_fq2,
    output:
        sorted_bam="output/wasp/remapped/{sampleName}_remapped_sorted.bam",
        bai="output/wasp/remapped/{sampleName}_remapped_sorted.bam.bai",
        stats="output/wasp/remapped/{sampleName}_stats.txt" 
    params:
        bwa_index=config['bwa_index'],
        bwa_version=config['bwaVers'],
        samtools_version=config['samtoolsVers']
    threads: 8
    log:
        err="output/logs/remap_wasp_{sampleName}.err",
        out="output/logs/remap_wasp_{sampleName}.out"
    shell:
        """
        module load bwa/{params.bwa_version}
        module load samtools/{params.samtools_version}
        mkdir -p output/wasp/remapped
        
        # Remap with BWA and pipe to sorting
        bwa mem -t {threads} -M {params.bwa_index} {input.fq1} {input.fq2} | \
        samtools view -b | \
        samtools sort -o {output.sorted_bam} \
        1> {log.out} 2> {log.err}

        samtools flagstat {output.sorted_bam} > {output.stats} 2>> {log.err}

        # Index the sorted BAM
        samtools index {output.sorted_bam} 1>> {log.out} 2>> {log.err}
        """

rule filter_remapped_reads:
    input:
        to_remap_bam=rules.find_intersecting_snps.output.to_remap_bam,
        remapped_bam=rules.remap_wasp_reads.output.sorted_bam,
        remapped_bai=rules.remap_wasp_reads.output.bai
    output:
        keep_bam="output/wasp/filtered_bam/{sampleName}.keep.bam"
    params:
        wasp_path=config['wasp_dir'],
        wasp_version=config['waspVer']
    log:
        "output/logs/filter_remapped_{sampleName}.log"
    shell:
        """
        module add wasp/{params.wasp_version}
        
        mkdir -p output/wasp/filtered_bam
        
        python {params.wasp_path}/mapping/filter_remapped_reads.py \
            {input.to_remap_bam} \
            {input.remapped_bam} \
            {output.keep_bam} \
            2> {log}
        """

rule merge:
    input:
        keepbam=rules.filter_remapped_reads.output.keep_bam,
        interbam=rules.find_intersecting_snps.output.keep_bam
    output:
        merged_bam="output/wasp/merge/{sampleName}.keep.merge.bam",
        sort_merged_bam="output/wasp/merge/{sampleName}.sort_keep.merge.bam",
        sort_merged_bai="output/wasp/merge/{sampleName}.sort_keep.merge.bam.bai"
    params:
        samtools_version=config['samtoolsVers']
    threads: 4
    log:
        "output/logs/merge_{sampleName}.log"
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p output/wasp/merge

        samtools merge {output.merged_bam} {input.keepbam} {input.interbam}
        samtools sort -o {output.sort_merged_bam} {output.merged_bam}
        samtools index {output.sort_merged_bam}

        """

rule rmDupWASP:
    input:
        merged_bam=rules.merge.output.sort_merged_bam
    output:
        rmdup_bam="output/wasp/rmdup/{sampleName}.rmdup.bam",
        sort_rmdup_bam="output/wasp/rmdup/{sampleName}.sorted_wasp.bam",
        index_bam="output/wasp/rmdup/{sampleName}.sorted_wasp.bam.bai",
        stats="output/wasp/rmdup/{sampleName}_stats.txt"
    params:
        pythonVer=config['python'],
        wasp_path=config['wasp_dir'],
        wasp_version=config['waspVer'],
        samtools_version=config['samtoolsVers']
    threads:8
    log:
        err="output/logs/rmdup_wasp_{sampleName}.err",
        out="output/logs/rmdup_wasp_{sampleName}.out"
    shell:
        """
        module add wasp/{params.wasp_version}
        module load samtools/{params.samtools_version}
        module load python/{params.pythonVer}
        mkdir -p output/wasp/rmdup

        python {params.wasp_path}/mapping/rmdup_pe.py {input.merged_bam} {output.rmdup_bam}

        samtools sort -o {output.sort_rmdup_bam} {output.rmdup_bam} 1> {log.out} 2> {log.err}
        samtools index {output.sort_rmdup_bam} 1>> {log.out} 2>> {log.err}
        samtools flagstat {output.sort_rmdup_bam} > {output.stats} 2>> {log.err}

        """

rule rmblacklist_region:
    input:
        bam=rules.rmDupWASP.output.sort_rmdup_bam
    output:
        final_bam="output/wasp/blk_filter/{sampleName}.blkfilter.bam",
        sort_bam="output/wasp/blk_filter/{sampleName}.sorted_final.bam",
        index_bam="output/wasp/blk_filter/{sampleName}.sorted_final.bam.bai",
        stats="output/wasp/blk_filter/{sampleName}_stats.txt"
    threads:4
    params:
        bedtools_ver=config['bedtools'],
        samtools_version=config['samtoolsVers'],
        Blacklist_bed=config['blacklist']
    log:
        err="output/logs/blklistFilter_wasp_{sampleName}.err",
        out="output/logs/blklistFilter_wasp_{sampleName}.out"
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        module load samtools/{params.samtools_version}
        mkdir -p output/wasp/blk_filter

        bedtools intersect -v -abam {input.bam} -b {params.Blacklist_bed} > {output.final_bam} 
        
        samtools sort -o {output.sort_bam} {output.final_bam}
        samtools index {output.sort_bam} 1>> {log.out} 2>> {log.err}
        samtools flagstat {output.sort_bam} > {output.stats} 2>> {log.err}

        """

rule run_ataqv:
    input:
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        txt="output/wasp/ataqv/{sampleName}.ataqv.txt",
        json="output/wasp/ataqv/{sampleName}.ataqv.json"
    params:
        tss_file="/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.tss.bed.gz",
        ataqv_version=config["ataqvVers"],
        tss_extension=config["tssextension"]
    threads: 6
    log:
        err="output/logs/ataqv_afterWASP_{sampleName}.err",
        out="output/logs/ataqv_afterWASP_{sampleName}.out"
    shell:
        """
        module load ataqv/{params.ataqv_version}
        mkdir -p output/wasp/ataqv

        ataqv \
            --tss-file {params.tss_file} \
            --tss-extension {params.tss_extension} \
            --metrics-file {output.json} \
            --ignore-read-groups human \
            {input.bam} > {output.txt} \
            2> {log.err}
         #mkarv /work/users/s/e/seyoun/CQTL_AI/output/ataqv/final_atacQC_summary *.json   
        """


rule bamQC:
    input:
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        qc="output/wasp/bamQC/{sampleName}_bamqQC"
    params:
        script="/work/users/s/e/seyoun/CQTL_AI/scripts/ATACseq_QC.R",
        r_version=config["r"]
    log:
        err="output/logs/atacseqQC_afterWASP_{sampleName}.err",
        out="output/logs/atacseqQC_afterWASP_{sampleName}.out"
    threads: 6
    shell:
        """
        module load r/{params.r_version}
        mkdir -p output/wasp/bamQC
        Rscript {params.script} {input.bam} output/wasp/bamQC/ .sorted_final.bam  1> {log.out} 2> {log.err}
        """

rule macs3_peakcalling:
    input:
        bam = rules.rmblacklist_region.output.sort_bam
    output:
        macs3="output/peaks/macs3/{sampleName}_peaks.narrowPeak"
    log:
        err="output/logs/macs3_{sampleName}.err"
    threads:2
    params:
        python_ver = config['python']
    shell:
        """
        module load python/{params.python_ver}
        source /users/s/e/seyoun/tools/MACS3env/bin/activate

        mkdir -p output/peaks/macs3/

        macs3 callpeak -t {input.bam} \
                -f BAMPE \
                -g hs \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                --keep-dup all \
                -q 0.01 \
                --call-summits \
                -n {wildcards.sampleName} \
                --outdir output/peaks/macs3/ 2> {log.err}
        """

rule hmmratac_peaks:
    input:
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        peak="output/peaks/hmmratac/{sampleName}_accessible_regions.narrowPeak"
    log:
        err='output/logs/hmmratac_{sampleName}.err'
    params:
        python_ver=config['python']
    threads:2
    shell:
        """
        module load python/{params.python_ver}
        source /users/s/e/seyoun/tools/MACS3env/bin/activate

        mkdir -p output/peaks/hmmratac/

         macs3 hmmratac \
            -i {input.bam} \
            -f BAMPE \
            --keep-duplicates \
            --minlen 100 \
            --save-states \
            --binsize 50 \
            -n {wildcards.sampleName} \
            --outdir output/peaks/hmmratac/ 2> {log.err}

        """
rule macs2_peaks:
    input:
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        peak="output/peaks/macs2/{sampleName}_peaks.narrowPeak"
    log:
        err='output/logs/macs2_{sampleName}.err'
    threads:2
    params:
        python_ver=config['python'],
        macs2_ver=config['macs2']
    shell:
        """
        module load python/{params.python_ver}
        module load macs/{params.macs2_ver}

        mkdir -p output/peaks/macs2

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
            --outdir output/peaks/macs2 2> {log.err}
        """
rule rocco_peaks:
    input:
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        peak="output/peaks/rocco/{sampleName}_peaks.bed"
    log:
        err='output/logs/rocco_{sampleName}.err'
    params:
        python_ver=config['python'],
        bedtools_ver=config['bedtools'],
        genomesize=config['effective_genome_size'],
        chrsize=config['chrsize']
    threads:4
    shell:
        """
        module load python/{params.python_ver}
        module load bedtools/{params.bedtools_ver}
        mkdir -p output/peaks/rocco

        rocco -i {input.bam} \
                --effective_genome_size {params.genomesize} \
                -o {output.peak} \
                --chrom_sizes_file {params.chrsize} \
                --step 10 \
                  --transform_logpc \
                  --transform_logpc_const 1.0 \
                  --transform_medfilt \
                  --transform_medfilt_kernel 150 \
                  --gamma 2.0 \
                  --budget 0.03 \
                  --min_mapping_score 30 \
                  --flag_include 66 \
                  --flag_exclude 1284 \
                  --threads 4 \
                  --verbose

        """

rule macs3_peak_count:
    input:
        macs3_peak = lambda wildcards: expand("output/peaks/macs3/{sp_nm_cond}_peaks.narrowPeak", sp_nm_cond=condition_samples[wildcards.cond]),
        bam=lambda wildcards: expand("output/wasp/blk_filter/{sp_nm_cond}.sorted_final.bam",sp_nm_cond=condition_samples[wildcards.cond])
    output:
        macs3_bed="output/peaks/merged_cond/{cond}_macs3_merged.bed",
        macs3_saf="output/peaks/merged_cond/{cond}_macs3_merged.saf",
        count="output/peaks/merged_cond/{cond}_macs3_merged_counts.txt"
    log:
        err="output/logs/{cond}_macs3_merge_peakcount.err",
        out="output/logs/{cond}_macs3_merge_peakcount.out"
    threads: 8
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread'],
        min_samples=6
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/merged_cond

        multiIntersectBed -i {input.macs3_peak} | \
        awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
        awk '($3-$2) >= 50' | \
        bedtools merge > {output.macs3_bed} 2> {log.err}

        #Convert to the saf

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.macs3_bed} > {output.macs3_saf}

         #Count reads in the region
         
         featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.macs3_saf} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}


        """

rule macs3_merge_conditions:
    input:
        bed=expand("output/peaks/merged_cond/{cond}_macs3_merged.bed", cond=['CTL', 'FNF']),
        bam=expand("output/wasp/blk_filter/{sampleName}.sorted_final.bam", sampleName=bam_files.keys())
    output:
        bed="output/peaks/merged/allsamples_macs3_merged.bed",
        saf="output/peaks/merged/allsamples_macs3_merged.saf",
        counts="output/peaks/merged/allsamples_macs3_merged_counts.txt"
    log:
        err="output/logs/macs3_merged_counts_Allsamples.err",
        out="output/logs/macs3_merged_counts_Allsamples.out"
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        
        mkdir -p output/peaks/merged

        cat {input.bed} |sort -k1,1 -k2,2n | bedtools merge > {output.bed} 2>> {log.err}

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf} 

        
        #Count reads in the region
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.saf} \
                -o {output.counts} \
                {input.bam} 2>> {log.err} 1>> {log.out}


        """

rule macs3_frip_calc:
    input:
        mergedSAF=rules.macs3_merge_conditions.output.saf,
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        count="output/peaks/macs3/indv/{sampleName}_count.txt"
    log:
        err="output/logs/macs3_{sampleName}_count.err",
        out="output/logs/macs3_{sampleName}_count.out"
    threads: 4
    params:
        subreadVer=config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/macs3/indv

        #count-indv
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {input.mergedSAF} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """

rule macs2_peak_count:
    input:
        peak = lambda wildcards: expand("output/peaks/macs2/{sp_nm_cond}_peaks.narrowPeak", sp_nm_cond=condition_samples[wildcards.cond]),
        bam = lambda wildcards: expand("output/wasp/blk_filter/{sp_nm_cond}.sorted_final.bam",sp_nm_cond=condition_samples[wildcards.cond])
    output:
        bed="output/peaks/merged_cond/{cond}_macs2_merged.bed",
        saf="output/peaks/merged_cond/{cond}_macs2_merged.saf",
        count="output/peaks/merged_cond/{cond}_macs2_merged_counts.txt"
    log:
        err="output/logs/{cond}_macs2_merge_peakcount.err",
        out="output/logs/{cond}_macs2_merge_peakcount.out"
    threads: 8
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread'],
        min_samples=6
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/merged_cond

        multiIntersectBed -i {input.peak} | \
        awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
        awk '($3-$2) >= 50' | \
        bedtools merge > {output.bed} 2> {log.err}

        #Convert to the saf

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}


        #Count reads in the region
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.saf} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        #--primary only count primary alignment
        #-C do not count reads where the pairs are mapped to different chromosomes

        """
rule macs2_merge_conditions:
    input:
        bed=expand("output/peaks/merged_cond/{cond}_macs2_merged.bed", cond=['CTL', 'FNF']),
        bam=expand("output/wasp/blk_filter/{sampleName}.sorted_final.bam", sampleName=bam_files.keys())
    output:
        bed="output/peaks/merged/allsamples_macs2_merged.bed",
        saf="output/peaks/merged/allsamples_macs2_merged.saf",
        counts="output/peaks/merged/allsamples_macs2_merged_counts.txt"
    log:
        err="output/logs/macs2_merged_counts_Allsamples.err",
        out="output/logs/macs2_merged_counts_Allsamples.out"
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        
        mkdir -p output/peaks/merged

        cat {input.bed} |sort -k1,1 -k2,2n | bedtools merge > {output.bed} 2>> {log.err}

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf} 

        
        #Count reads in the region
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.saf} \
                -o {output.counts} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """

rule macs2_frip_calc:
    input:
        mergedSAF=rules.macs2_merge_conditions.output.saf,
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        count="output/peaks/macs2/indv/{sampleName}_count.txt"
    log:
        err="output/logs/macs2_{sampleName}_count.err",
        out="output/logs/macs2_{sampleName}_count.out"
    threads: 4
    params:
        subreadVer=config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/macs2/indv

        #count-indv
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {input.mergedSAF} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """

rule hmmratac_peak_count:
    input:
        peak = lambda wildcards: expand("output/peaks/hmmratac/{sp_nm_cond}_accessible_regions.narrowPeak", sp_nm_cond=condition_samples[wildcards.cond]),
        bam = lambda wildcards: expand("output/wasp/blk_filter/{sp_nm_cond}.sorted_final.bam",sp_nm_cond=condition_samples[wildcards.cond])
    output:
        bed="output/peaks/merged_cond/{cond}_hmmratac_merged.bed",
        saf="output/peaks/merged_cond/{cond}_hmmratac_merged.saf",
        count="output/peaks/merged_cond/{cond}_hmmratac_merged_counts.txt"
    log:
        err="output/logs/{cond}_hmmratac_merge_peakcount.err",
        out="output/logs/{cond}_hmmratac_merge_peakcount.out"
    threads: 8
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread'],
        min_samples=6
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/merged_cond

        multiIntersectBed -i {input.peak} | \
        awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
        awk '($3-$2) >= 40' | \
        bedtools merge > {output.bed} 2> {log.err}

        #Convert to the saf

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}


        #Count reads in the region
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.saf} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        #--primary only count primary alignment
        #-C do not count reads where the pairs are mapped to different chromosomes

        """

rule hmmratac_merge_conditions:
    input:
        bed=expand("output/peaks/merged_cond/{cond}_hmmratac_merged.bed", cond=['CTL', 'FNF']),
        bam=expand("output/wasp/blk_filter/{sampleName}.sorted_final.bam", sampleName=bam_files.keys())
    output:
        bed="output/peaks/merged/allsamples_hmmratac_merged.bed",
        saf="output/peaks/merged/allsamples_hmmratac_merged.saf",
        counts="output/peaks/merged/allsamples_hmmratac_merged_counts.txt"
    log:
        err="output/logs/hmmratac_merged_counts_Allsamples.err",
        out="output/logs/hmmratac_merged_counts_Allsamples.out"
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}

        mkdir -p output/peaks/merged

        cat {input.bed} |sort -k1,1 -k2,2n | bedtools merge > {output.bed} 2>> {log.err}

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}


        #Count reads in the region
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.saf} \
                -o {output.counts} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """

rule hmmratac_frip_calc:
    input:
        mergedSAF=rules.hmmratac_merge_conditions.output.saf,
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        count="output/peaks/hmmratac/indv/{sampleName}_count.txt"
    log:
        err="output/logs/hmmratac_{sampleName}_count.err",
        out="output/logs/hmmratac_{sampleName}_count.out"
    threads: 4
    params:
        subreadVer=config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/hmmratac/indv

        #count-indv
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {input.mergedSAF} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """
rule rocco_peak_count:
    input:
        peak = lambda wildcards: expand("output/peaks/rocco/{sp_nm_cond}_peaks.bed", sp_nm_cond=condition_samples[wildcards.cond]),
        bam = lambda wildcards: expand("output/wasp/blk_filter/{sp_nm_cond}.sorted_final.bam",sp_nm_cond=condition_samples[wildcards.cond])
    output:
        bed="output/peaks/merged_cond/{cond}_rocco_merged.bed",
        saf="output/peaks/merged_cond/{cond}_rocco_merged.saf",
        count="output/peaks/merged_cond/{cond}_rocco_merged_counts.txt"
    log:
        err="output/logs/{cond}_rocco_merge_peakcount.err",
        out="output/logs/{cond}_rocco_merge_peakcount.out"
    threads: 8
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread'],
        min_samples=6
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/merged_cond

        multiIntersectBed -i {input.peak} | \
        awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
        awk '($3-$2) >= 40' | \
        bedtools merge > {output.bed} 2> {log.err}

        #Convert to the saf

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}


        #Count reads in the region
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.saf} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        #--primary only count primary alignment
        #-C do not count reads where the pairs are mapped to different chromosomes
        """

rule rocco_merge_conditions:
    input:
        bed=expand("output/peaks/merged_cond/{cond}_rocco_merged.bed", cond=['CTL', 'FNF']),
        bam=expand("output/wasp/blk_filter/{sampleName}.sorted_final.bam", sampleName=bam_files.keys())
    output:
        bed="output/peaks/merged/allsamples_rocco_merged.bed",
        saf="output/peaks/merged/allsamples_rocco_merged.saf",
        counts="output/peaks/merged/allsamples_rocco_merged_counts.txt"
    log:
        err="output/logs/rocco_merged_counts_Allsamples.err",
        out="output/logs/rocco_merged_counts_Allsamples.out"
    params:
        bedtoolsVer=config['bedtools'],
        subreadVer=config['subread']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}

        mkdir -p output/peaks/merged

        cat {input.bed} |sort -k1,1 -k2,2n | bedtools merge > {output.bed} 2>> {log.err}

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
         {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}


        #Count reads in the region
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {output.saf} \
                -o {output.counts} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """
rule rocco_frip_calc:
    input:
        mergedSAF=rules.rocco_merge_conditions.output.saf,
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        count="output/peaks/rocco/indv/{sampleName}_count.txt"
    log:
        err="output/logs/rocco_{sampleName}_count.err",
        out="output/logs/rocco_{sampleName}_count.out"
    threads: 4
    params:
        subreadVer=config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/rocco/indv

        #count-indv
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {input.mergedSAF} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """

rule signal:
    input:
        bam = rules.rmblacklist_region.output.sort_bam
    output:
        signal = "output/signals/indv/{sampleName}.bw"
    log:
        err = 'output/logs/signal_{sampleName}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        binSize=config['bin_size'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    threads:4
    shell:
        """ 
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/indv/

        bamCoverage \
            --bam {input.bam} \
            -o {output.signal} \
            --binSize {params.binSize} \
            --normalizeUsing {params.NormOption} \
            --effectiveGenomeSize {params.effective_genomeSize} \
            --extendReads \
            --minMappingQuality 30 \
            --samFlagInclude 66 \
            --samFlagExclude 1284 \
            --ignoreForNormalization chrX chrY chrM \
            --numberOfProcessors 4
            > {log.err} 2>&1

        """

# First define the empty dictionaries
ctl_bams = {}
fnf_bams = {}
for idx, sample in enumerate(samples['mn']):
    bam_path= f"output/wasp/blk_filter/{sample}.sorted_final.bam"
    if 'CTL' in sample:
        ctl_bams[sample] = bam_path
    if 'FNF' in sample:
        fnf_bams[sample] = bam_path

print("CTL BAMs:", ctl_bams)
print("FNF BAMs:", fnf_bams)

condition_bams = {
   'CTL': list(ctl_bams.values()),
   'FNF': list(fnf_bams.values())
}

rule merge_by_condition:
    input:
        bams = lambda wildcards: condition_bams[wildcards.condition]
    output:
        merged_bam = "output/wasp/merged/{condition}_merged.bam",
        merged_bai = "output/wasp/merged/{condition}_merged.bam.bai",
        stats = "output/wasp/merged/{condition}_merged_stats.txt"
    params:
        samtools_version = config['samtoolsVers']
    threads: 8
    log:
        err = "output/logs/merge_{condition}.err",
        out = "output/logs/merge_{condition}.out"
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p output/wasp/merged

        # Merge BAM files
        samtools merge -@ {threads} {output.merged_bam} {input.bams} 1>> {log.out} 2>> {log.err}

        # Index merged BAM
        samtools index {output.merged_bam} 1>> {log.out} 2>> {log.err}

        # Get stats for merged file
        samtools flagstat {output.merged_bam} > {output.stats} 2>> {log.err}
        """

rule mergeSignal:
    input:
        bam = rules.merge_by_condition.output.merged_bam
    output:
        signal='output/signals/merged_signal/{condition}_merged.bw'
    log:
        err='output/logs/mergedSignal_{condition}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        binSize=config['bin_size'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/merged_signal

        bamCoverage \
            --bam {input.bam} \
            -o {output.signal} \
            --binSize {params.binSize} \
            --normalizeUsing {params.NormOption} \
            --effectiveGenomeSize {params.effective_genomeSize} \
            --extendReads \
            --minMappingQuality 30 \
            --samFlagInclude 66 \
            --samFlagExclude 1284 \
            --ignoreForNormalization chrX chrY chrM \
            --numberOfProcessors 4
            > {log.err} 2>&1
        """



