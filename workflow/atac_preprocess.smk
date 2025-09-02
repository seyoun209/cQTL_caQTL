#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

# -- Load samplesheet and build file lists --
S = pd.read_csv(config["samplesheet"], sep="\t", dtype=str).fillna("")

# Guess Read1/Read2 if missing (optional)
if "Read1" not in S.columns or "Read2" not in S.columns:
    def guess(r):
        last = r.iloc[-1]
        if isinstance(last, str) and last.endswith(".fastq.gz"):
            r1 = os.path.join(r.get("Sequencing_Directory",""), last) if not last.startswith("/") else last
            r2 = r1.replace("_R1_", "_R2_")
            return pd.Series({"Read1": r1, "Read2": r2})
        return pd.Series({"Read1":"", "Read2":""})
    S[["Read1","Read2"]] = S.apply(guess, axis=1)

# Make paths absolute using Sequencing_Directory if needed
for col in ("Read1","Read2"):
    S[col] = S.apply(
        lambda r: r[col] if not r.get("Sequencing_Directory") or r[col].startswith("/")
        else os.path.join(r["Sequencing_Directory"], r[col]),
        axis=1
    )

merge_keys = config.get("mergeBy", config.get("fileNamesFrom", ["Proj","Donor","Condition","Tissue","Protocol_notes"]))
S["mn"] = S[merge_keys].agg("_".join, axis=1).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True)

READ1 = S.groupby("mn")["Read1"].apply(list).to_dict()
READ2 = S.groupby("mn")["Read2"].apply(list).to_dict()
SAMPLES = sorted(READ1.keys())

onsuccess:
    print("ATAC_pre_process completed successfully! Wahoo!")

rule all:
    input:
        # fastqs
        [OUT("fastq", f"{s}_R1.fastq.gz") for s in SAMPLES],
        [OUT("fastq", f"{s}_R2.fastq.gz") for s in SAMPLES],
        # fastqc
        [OUT("QC", f"{s}_R1_fastqc.zip") for s in SAMPLES],
        [OUT("QC", f"{s}_R1_fastqc.html") for s in SAMPLES],
        [OUT("QC", f"{s}_R2_fastqc.zip") for s in SAMPLES],
        [OUT("QC", f"{s}_R2_fastqc.html") for s in SAMPLES],
        # trim
        [OUT("trim", f"{s}_R1_val_1.fq.gz") for s in SAMPLES],
        [OUT("trim", f"{s}_R2_val_2.fq.gz") for s in SAMPLES],
        [OUT("trim", f"{s}_R1.fastq.gz_trimming_report.txt") for s in SAMPLES],
        [OUT("trim", f"{s}_R2.fastq.gz_trimming_report.txt") for s in SAMPLES],
        # align/dedup/filter
        [OUT("align", f"{s}_sorted.bam") for s in SAMPLES],
        [OUT("filtered/blk_filter", f"{s}.sorted_final.bam") for s in SAMPLES],
        # metrics
        [OUT("metrics", f"{s}_insert_size_metrics.txt") for s in SAMPLES],
        [OUT("metrics", f"{s}_insert_size_histogram.pdf") for s in SAMPLES],
        # QC
        [OUT("ataqv", f"{s}.ataqv.json") for s in SAMPLES],
        [OUT("ataqv", f"{s}.ataqv.txt") for s in SAMPLES],
        [OUT("bamQC", f"{s}_bamQC") for s in SAMPLES],
        # bigwigs
        [OUT("signals/indv", f"{s}.bw") for s in SAMPLES],
        # bam samplesheet
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
        errR1=OUT("logs", "{sample}_R1_catReads.err"),
        errR2=OUT("logs", "{sample}_R2_catReads.err")
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log.errR1}
        cat {input.R2} > {output.R2} 2> {log.errR2}
        """

rule fastqc:
    input:
        R1=rules.catReads.output.R1,
        R2=rules.catReads.output.R2
    output:
        zip1=OUT("QC", "{sample}_R1_fastqc.zip"),
        html1=OUT("QC", "{sample}_R1_fastqc.html"),
        zip2=OUT("QC", "{sample}_R2_fastqc.zip"),
        html2=OUT("QC", "{sample}_R2_fastqc.html")
    threads: 2
    params:
        version=config["fastqcVers"]
    log:
        err=OUT("logs", "fastqc_{sample}.err"),
        out=OUT("logs", "fastqc_{sample}.out")
    shell:
        """
        module load fastqc/{params.version}
        mkdir -p {OUT("QC")}
        fastqc {input.R1} {input.R2} -o {OUT("QC")} 1> {log.out} 2> {log.err}
        """

rule trim:
    input:
        R1=rules.catReads.output.R1,
        R2=rules.catReads.output.R2
    output:
        trim1=OUT("trim", "{sample}_R1_val_1.fq.gz"),
        trim2=OUT("trim", "{sample}_R2_val_2.fq.gz"),
        report1=OUT("trim", "{sample}_R1.fastq.gz_trimming_report.txt"),
        report2=OUT("trim", "{sample}_R2.fastq.gz_trimming_report.txt")
    threads: 4
    params:
        version=config["trim_galore"]
    log:
        err=OUT("logs", "trim_{sample}.err"),
        out=OUT("logs", "trim_{sample}.out")
    shell:
        """
        module load trim_galore/{params.version}
        module load pigz
        mkdir -p {OUT("trim")}
        trim_galore -o {OUT("trim")} --cores {threads} --paired {input.R1} {input.R2} 2> {log.err}
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
        samtools_version=config["samtoolsVers"]
    threads: 8
    log:
        err=OUT("logs", "align_{sample}.err"),
        out=OUT("logs", "align_{sample}.out")
    shell:
        """
        module load bwa/{params.bwa_version}
        module load samtools/{params.samtools_version}
        mkdir -p {OUT("align")}
        bwa mem -t {threads} -M {params.bwa_index} {input.trim1} {input.trim2} | \
          samtools view -q 30 -b | \
          samtools sort -o {output.sortedBam} 1> {log.out} 2> {log.err}
        samtools flagstat {output.sortedBam} > {output.stats} 2>> {log.err}
        """

rule mark_duplicates:
    input:
        bam=rules.align.output.sortedBam
    output:
        dedup_bam=OUT("dedup", "{sample}_dedup.bam"),
        metrics=OUT("dedup", "{sample}_dup_metrics.txt"),
        index=OUT("dedup", "{sample}_dedup.bam.bai")
    threads: 4
    params:
        picard_version=config["picardVers"],
        java_version=config["javaVers"],
        samtools_version=config["samtoolsVers"]
    log:
        err=OUT("logs", "mark_duplicates_{sample}.err"),
        out=OUT("logs", "mark_duplicates_{sample}.out")
    shell:
        """
        module load java/{params.java_version}
        module load picard/{params.picard_version}
        module load samtools/{params.samtools_version}
        mkdir -p {OUT("dedup")}
        java -Xmx16g -jar /nas/longleaf/apps/picard/{params.picard_version}/picard-{params.picard_version}/picard.jar MarkDuplicates \
            I={input.bam} O={output.dedup_bam} M={output.metrics} REMOVE_DUPLICATES=true \
            1> {log.out} 2> {log.err}
        samtools index {output.dedup_bam} 1>> {log.out} 2>> {log.err}
        """

rule filter_mitochondrial_reads:
    input:
        dedup_bam=rules.mark_duplicates.output.dedup_bam,
        index=rules.mark_duplicates.output.index
    output:
        filtered_bam=OUT("filtered", "{sample}_filtered.bam"),
        index=OUT("filtered", "{sample}_filtered.bam.bai")
    threads: 4
    params:
        samtools_version=config["samtoolsVers"]
    log:
        err=OUT("logs", "filter_mitochondrial_{sample}.err"),
        out=OUT("logs", "filter_mitochondrial_{sample}.out")
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p {OUT("filtered")}
        samtools view -bh {input.dedup_bam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
          -o {output.filtered_bam} 1> {log.out} 2> {log.err}
        samtools index {output.filtered_bam} 1>> {log.out} 2>> {log.err}
        """

rule rmblacklist_region:
    input:
        bam=rules.filter_mitochondrial_reads.output.filtered_bam
    output:
        final_bam=OUT("filtered/blk_filter", "{sample}.blkfilter.bam"),
        sort_bam=OUT("filtered/blk_filter", "{sample}.sorted_final.bam"),
        index_bam=OUT("filtered/blk_filter", "{sample}.sorted_final.bam.bai"),
        stats=OUT("filtered/blk_filter", "{sample}_stats.txt")
    threads: 4
    params:
        bedtools_ver=config["bedtools"],
        samtools_version=config["samtoolsVers"],
        Blacklist_bed=config["blacklist"]
    log:
        err=OUT("logs", "blk_filter_{sample}.err"),
        out=OUT("logs", "blk_filter_{sample}.out")
    shell:
        """
        module load bedtools/{params.bedtools_ver}
        module load samtools/{params.samtools_version}
        mkdir -p {OUT("filtered/blk_filter")}
        bedtools intersect -v -abam {input.bam} -b {params.Blacklist_bed} > {output.final_bam}
        samtools sort -o {output.sort_bam} {output.final_bam}
        samtools index {output.sort_bam} 1>> {log.out} 2>> {log.err}
        samtools flagstat {output.sort_bam} > {output.stats} 2>> {log.err}
        """

rule collect_insert_size_metrics:
    input:
        bam=rules.rmblacklist_region.output.final_bam
    output:
        metrics=OUT("metrics", "{sample}_insert_size_metrics.txt"),
        histogram=OUT("metrics", "{sample}_insert_size_histogram.pdf")
    threads: 4
    params:
        picard_version=config["picardVers"],
        java_version=config["javaVers"],
        R_version=config["r"]
    log:
        err=OUT("logs", "insert_size_metrics_{sample}.err"),
        out=OUT("logs", "insert_size_metrics_{sample}.out")
    shell:
        """
        module load java/{params.java_version}
        module load picard/{params.picard_version}
        module load r/{params.R_version}
        mkdir -p {OUT("metrics")}
        java -jar /nas/longleaf/apps/picard/{params.picard_version}/picard-{params.picard_version}/picard.jar CollectInsertSizeMetrics \
          I={input.bam} O={output.metrics} H={output.histogram} M=0.05 ASSUME_SORTED=true \
          1> {log.out} 2> {log.err}
        """

rule add_read_groups:
    input:
        bam=rules.rmblacklist_region.output.final_bam,
        bai=rules.rmblacklist_region.output.index_bam
    output:
        bam=OUT("align", "{sample}_RG.bam"),
        bai=OUT("align", "{sample}_RG.bam.bai")
    params:
        picard_version=config["picardVers"],
        samtools_version=config["samtoolsVers"]
    threads: 4
    log:
        err=OUT("logs", "add_read_groups_{sample}.err"),
        out=OUT("logs", "add_read_groups_{sample}.out")
    shell:
        """
        module load picard/{params.picard_version}
        module load samtools/{params.samtools_version}
        mkdir -p {OUT("align")}
        picard AddOrReplaceReadGroups I={input.bam} O={output.bam} RGSM={wildcards.sample} RGPL=ILLUMINA RGLB=lib1 RGPU=unit1 \
          1> {log.out} 2> {log.err}
        samtools index {output.bam} 1>> {log.out} 2>> {log.err}
        """

rule run_ataqv:
    input:
        bam=rules.add_read_groups.output.bam
    output:
        txt=OUT("ataqv", "{sample}.ataqv.txt"),
        json=OUT("ataqv", "{sample}.ataqv.json")
    params:
        tss_file="/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.tss.bed.gz",
        ataqv_version=config["ataqvVers"],
        tss_extension=config["tssextension"]
    threads: 6
    log:
        err=OUT("logs", "ataqv_{sample}.err"),
        out=OUT("logs", "ataqv_{sample}.out")
    shell:
        """
        module load ataqv/{params.ataqv_version}
        mkdir -p {OUT("ataqv")}
        ataqv --tss-file {params.tss_file} --tss-extension {params.tss_extension} \
          --metrics-file {output.json} --ignore-read-groups human {input.bam} > {output.txt} 2> {log.err}
        """

rule qc_r:
    input:
        bam=rules.add_read_groups.output.bam
    output:
        qc=OUT("bamQC", "{sample}_bamQC")
    params:
        script="/work/users/s/e/seyoun/CQTL_AI/scripts/ATACseq_QC.R",
        r_version=config["r"]
    threads: 6
    log:
        err=OUT("logs", "atacseqQC_{sample}.err"),
        out=OUT("logs", "atacseqQC_{sample}.out")
    shell:
        """
        module load r/{params.r_version}
        Rscript {params.script} {input.bam} "{OUT('bamQC')}" "_RG.bam" 1> {log.out} 2> {log.err}
        """

rule generate_bam_samplesheet:
    input:
        bam=[OUT("filtered/blk_filter", f"{s}.sorted_final.bam") for s in SAMPLES]
    output:
        bam_samplesheet=OUT("bam_atac_CQTL_all_samplesheet.txt")
    run:
        df = pd.read_csv(config["samplesheet"], sep="\t")
        df["Bam_file"] = df.apply(lambda r: f"{r['Proj']}_{r['Donor']}_{r['Condition']}_{r['Tissue']}_{r['Protocol_notes']}.sorted_final.bam", axis=1)
        df["Bam_directory"] = OUT("filtered/blk_filter")
        df.to_csv(output.bam_samplesheet, sep="\t", index=False)

rule signal:
    input:
        bam=rules.rmblacklist_region.output.sort_bam
    output:
        signal=OUT("signals/indv", "{sample}.bw")
    threads: 4
    params:
        deeptools_ver=config["deeptoolsVers"],
        binSize=config["bin_size"],
        effective_genomeSize=config["effective_genome_size"],
        NormOption=config["normalize_option"]
    log:
        err=OUT("logs", "signal_{sample}.err")
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p {OUT("signals/indv")}
        bamCoverage --bam {input.bam} -o {output.signal} --binSize {params.binSize} \
          --normalizeUsing {params.NormOption} --effectiveGenomeSize {params.effective_genomeSize} \
          --extendReads --minMappingQuality 30 --samFlagInclude 66 --samFlagExclude 1284 \
          --ignoreForNormalization chrX chrY chrM --numberOfProcessors {threads} > {log.err} 2>&1
        """
