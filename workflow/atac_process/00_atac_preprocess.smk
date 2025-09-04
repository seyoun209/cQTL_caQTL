#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

# Load samplesheet and build file lists
S = pd.read_csv(config["samplesheet"], sep="\t", dtype=str).fillna("")

# Handle Read1/Read2 paths
if "Read1" not in S.columns or "Read2" not in S.columns:
    def guess(r):
        last = r.iloc[-1]
        if isinstance(last, str) and last.endswith(".fastq.gz"):
            r1 = os.path.join(r.get("Sequencing_Directory",""), last) if not last.startswith("/") else last
            r2 = r1.replace("_R1_", "_R2_")
            return pd.Series({"Read1": r1, "Read2": r2})
        return pd.Series({"Read1":"", "Read2":""})
    S[["Read1","Read2"]] = S.apply(guess, axis=1)

# Make paths absolute
for col in ("Read1","Read2"):
    S[col] = S.apply(
        lambda r: r[col] if not r.get("Sequencing_Directory") or r[col].startswith("/")
        else os.path.join(r["Sequencing_Directory"], r[col]),
        axis=1
    )

# Create sample names
merge_keys = config.get("fileNamesFrom", ["Proj","Donor","Condition","Tissue","Protocol_notes"])
S["mn"] = S[merge_keys].agg("_".join, axis=1).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True)

READ1 = S.groupby("mn")["Read1"].apply(list).to_dict()
READ2 = S.groupby("mn")["Read2"].apply(list).to_dict()
SAMPLES = sorted(READ1.keys())

onsuccess:
    print("ATAC preprocessing completed successfully!")
    print("Final BAM files: atac_output/filtered/blk_filter/*.sorted_final.bam")
    print("Ready for peak calling, signal tracks, or WASP processing")

rule all:
    input:
        # Main outputs - final processed BAMs
        [OUT("filtered/blk_filter", f"{s}.sorted_final.bam") for s in SAMPLES],
        [OUT("filtered/blk_filter", f"{s}.sorted_final.bam.bai") for s in SAMPLES],
        # QC outputs
        [OUT("ataqv", f"{s}.ataqv.json") for s in SAMPLES],
        [OUT("bamQC", f"{s}_bamQC") for s in SAMPLES],
        # Metrics
        [OUT("metrics", f"{s}_insert_size_metrics.txt") for s in SAMPLES],
        # Sample sheet for downstream processing
        OUT("bam_atac_CQTL_all_samplesheet.txt")

rule catReads:
    input:
        R1=lambda wc: READ1[wc.sample],
        R2=lambda wc: READ2[wc.sample]
    output:
        R1=OUT("fastq", "{sample}_R1.fastq.gz"),
        R2=OUT("fastq", "{sample}_R2.fastq.gz")
    threads: 2
    log:
        OUT("logs", "{sample}_catReads.log")
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log}
        cat {input.R2} > {output.R2} 2>> {log}
        """

rule trim:
    input:
        R1=rules.catReads.output.R1,
        R2=rules.catReads.output.R2
    output:
        trim1=OUT("trim", "{sample}_R1_val_1.fq.gz"),
        trim2=OUT("trim", "{sample}_R2_val_2.fq.gz")
    threads: 4
    params:
        version=config["trim_galore"],
        trim_dir=OUT('trim')
    log:
        OUT("logs", "trim_{sample}.log")
    shell:
        """
        module load trim_galore/{params.version}
        module load pigz
        mkdir -p {params.trim_dir}
        trim_galore -o {params.trim_dir} --cores {threads} --paired {input.R1} {input.R2} 2> {log}
        """

rule align:
    input:
        trim1=rules.trim.output.trim1,
        trim2=rules.trim.output.trim2
    output:
        sortedBam=OUT("align", "{sample}_sorted.bam"),
        stats=OUT("align", "{sample}_stats.txt")
    params:
        bwa_index=config["bwa_index"],
        bwa_version=config["bwaVers"],
        samtools_version=config["samtoolsVers"],
        align_dir=OUT("align")
    threads: 8
    log:
        OUT("logs", "align_{sample}.log")
    shell:
        """
        module load bwa/{params.bwa_version}
        module load samtools/{params.samtools_version}
        mkdir -p {params.align_dir}
        
        bwa mem -t {threads} -M {params.bwa_index} {input.trim1} {input.trim2} | \
          samtools view -q 30 -b | \
          samtools sort -o {output.sortedBam} 2> {log}
        samtools flagstat {output.sortedBam} > {output.stats} 2>> {log}
        """

rule add_read_groups:
    input:
        bam=rules.align.output.sortedBam
    output:
        bam=OUT("align", "{sample}_RG.bam"),
        bai=OUT("align", "{sample}_RG.bam.bai")
    params:
        picard_version=config["picardVers"],
        samtools_version=config["samtoolsVers"],
        java_version=config["javaVers"],
        align_dir=OUT("align")
    threads: 4
    log:
        OUT("logs", "add_read_groups_{sample}.log")
    shell:
        """
        module load java/{params.java_version}
        module load picard/{params.picard_version}
        module load samtools/{params.samtools_version}
        mkdir -p {params.align_dir}

        picard AddOrReplaceReadGroups I={input.bam} O={output.bam} \
          RGSM={wildcards.sample} RGPL=ILLUMINA RGLB=lib1 RGPU=unit1 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule mark_duplicates:
    input:
        bam=rules.add_read_groups.output.bam
    output:
        dedup_bam=OUT("dedup", "{sample}_dedup.bam"),
        metrics=OUT("dedup", "{sample}_dup_metrics.txt"),
        index=OUT("dedup", "{sample}_dedup.bam.bai")
    threads: 4
    params:
        picard_version=config["picardVers"],
        java_version=config["javaVers"],
        samtools_version=config["samtoolsVers"],
        dedup_dir=OUT("dedup")
    log:
        OUT("logs", "mark_duplicates_{sample}.log")
    shell:
        """
        module load java/{params.java_version}
        module load picard/{params.picard_version}
        module load samtools/{params.samtools_version}
        mkdir -p {params.dedup_dir}
        
        picard MarkDuplicates I={input.bam} O={output.dedup_bam} M={output.metrics} REMOVE_DUPLICATES=true 2> {log}
        samtools index {output.dedup_bam} 2>> {log}
        """


rule filter_mitochondrial_reads:
    input:
        dedup_bam=rules.mark_duplicates.output.dedup_bam
    output:
        filtered_bam=OUT("filtered", "{sample}_filtered.bam"),
        index=OUT("filtered", "{sample}_filtered.bam.bai")
    threads: 4
    params:
        samtools_version=config["samtoolsVers"],
        filtered_dir=OUT("filtered")
    log:
        OUT("logs", "filter_mitochondrial_{sample}.log")
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p {params.filtered_dir}
        
        samtools view -bh {input.dedup_bam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
          chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
          -o {output.filtered_bam} 2> {log}
        samtools index {output.filtered_bam} 2>> {log}
        """

rule remove_blacklist_regions:
    input:
        bam=rules.filter_mitochondrial_reads.output.filtered_bam
    output:
        final_bam=OUT("filtered/blk_filter", "{sample}.sorted_final.bam"),
        index_bam=OUT("filtered/blk_filter", "{sample}.sorted_final.bam.bai"),
        stats=OUT("filtered/blk_filter", "{sample}_stats.txt")
    threads: 4
    params:
        bedtools_ver=config["bedtools"],
        samtools_version=config["samtoolsVers"],
        blacklist_bed=config["blacklist"],
        blk_filter_dir=OUT("filtered/blk_filter")
    log:
        OUT("logs", "remove_blacklist_{sample}.log")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        module load samtools/{params.samtools_version}
        mkdir -p {params.blk_filter_dir}
        
        bedtools intersect -v -abam {input.bam} -b {params.blacklist_bed} | \
          samtools sort -o {output.final_bam} 2> {log}
        samtools index {output.final_bam} 2>> {log}
        samtools flagstat {output.final_bam} > {output.stats} 2>> {log}
        """

rule collect_insert_size_metrics:
    input:
        bam=rules.remove_blacklist_regions.output.final_bam
    output:
        metrics=OUT("metrics", "{sample}_insert_size_metrics.txt"),
        histogram=OUT("metrics", "{sample}_insert_size_histogram.pdf")
    threads: 4
    params:
        picard_version=config["picardVers"],
        java_version=config["javaVers"],
        R_version=config["r"],
        metrics_dir=OUT("metrics")
    log:
        OUT("logs", "insert_size_metrics_{sample}.log")
    shell:
        """
        module load java/{params.java_version}
        module load picard/{params.picard_version}
        module load r/{params.R_version}
        mkdir -p {params.metrics_dir}
        
        java -jar /nas/longleaf/apps/picard/{params.picard_version}/picard-{params.picard_version}/picard.jar CollectInsertSizeMetrics \
          I={input.bam} O={output.metrics} H={output.histogram} M=0.05 ASSUME_SORTED=true 2> {log}
        """

rule run_ataqv:
    input:
        bam=rules.remove_blacklist_regions.output.final_bam
    output:
        txt=OUT("ataqv", "{sample}.ataqv.txt"),
        json=OUT("ataqv", "{sample}.ataqv.json")
    params:
        tss_file=config["tss_annotation"],
        ataqv_bin=config["ataqvVers"],
        tss_extension=config["tssextension"],
        ataqv_dir=OUT("ataqv")
    threads: 6
    log:
        OUT("logs", "ataqv_{sample}.log")
    shell:
        """
        mkdir -p {params.ataqv_dir}
        
        {params.ataqv_bin} --tss-file {params.tss_file} --tss-extension {params.tss_extension} \
          --metrics-file {output.json} --ignore-read-groups human \
          {input.bam} > {output.txt} 2> {log}
        """

rule bam_qc:
    input:
        bam=rules.remove_blacklist_regions.output.final_bam
    output:
        qc=OUT("bamQC", "{sample}_bamQC")
    params:
        script=config["atac_qc_Rscript"],
        r_version=config["r"],
        bamqc_dir=OUT("bamQC")
    threads: 6
    log:
        OUT("logs", "bam_qc_{sample}.log")
    shell:
        """
        module load r/{params.r_version}
        mkdir -p {params.bamqc_dir}
        Rscript {params.script} {input.bam} "{params.bamqc_dir}" ".sorted_final.bam" 2> {log}
        """

rule generate_bam_samplesheet:
    input:
        bam=[OUT("filtered/blk_filter", f"{s}.sorted_final.bam") for s in SAMPLES]
    output:
        bam_samplesheet=OUT("bam_atac_CQTL_all_samplesheet.txt")
    run:
        df = pd.read_csv(config["samplesheet"], sep="\t")
        merge_keys = config.get("fileNamesFrom", ["Proj","Donor","Condition","Tissue","Protocol_notes"])
        df["mn"] = df[merge_keys].agg("_".join, axis=1).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True)
        df["Bam_file"] = df["mn"] + ".sorted_final.bam"
        df["Bam_directory"] = OUT("filtered/blk_filter")
        df.to_csv(output.bam_samplesheet, sep="\t", index=False)
