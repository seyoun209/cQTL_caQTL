#!/usr/bin/env python3

include: '/work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk'

# Define the list of all possible sample names
all_samples = [
    "CQTL_AM7778FE_FNF_Femur_1",
    "CQTL_AM7754_FNF_Ankle_replicate",
    "CQTL_AM7755_FNF_Ankle_replicate",
    "CQTL_AM7778FE_CTL_Femur_1",
    "CQTL_AM7754_CTL_Ankle_replicate",
    "CQTL_AM7755_CTL_Ankle_replicate"
    #"OTHER_SAMPLE_1",
    #"OTHER_SAMPLE_2"
]

# Filter samples to include only those with 'replicate' or 'Femur' in their names
selected_samples = [s for s in all_samples if "replicate" in s or "Femur" in s]

rule all_subset:
    input:
        "output/geno/pbs_geno_other/pbs.rename_wCHR.vcf.gz",
        "output/geno/fnf_geno_other/fnf.rename_wCHR.vcf.gz",
        "output/geno/pbs_geno_other/pbs.rename.vcf.gz",
        "output/geno/fnf_geno_other/fnf.rename.vcf.gz",
        expand("output/QC/{sampleName}_verifyBamID.{ext}",
            sampleName=selected_samples,
            ext=['selfSM','selfRG','bestRG','bestSM','depthRG','depthSM']),
        expand("output/geno/pbs_geno_other/wasp_vcf/chr{chrom}.pbs.rename.vcf.gz",
            chrom=range(1,23)),
        expand("output/geno/fnf_geno_other/wasp_vcf/chr{chrom}.fnf.rename.vcf.gz",
            chrom=range(1,23)),
        ataqv=expand("output/wasp/ataqv/{sampleName}.ataqv.json",
            sampleName=selected_samples),
        ataqv_txt=expand("output/wasp/ataqv/{sampleName}.ataqv.txt",
            sampleName=selected_samples),
        bamqc= expand("output/wasp/bamQC/{sampleName}_bamqQC", sampleName=selected_samples),
        # Individual files
        indv_counts = expand("output/peaks/{tool}/indv/{sampleName}_count.txt", tool=['macs2','macs3','hmmratac','rocco'], sampleName=selected_samples),
        indv_signal = expand("output/signals/indv/{sampleName}.bw", sampleName=selected_samples)


rule categorize_femur_replicate_samples:
    output:
        ctl_femur_replicate="output/geno/pbs_geno_other/CTL_samples_with_replacements.txt",
        fnf_femur_replicate="output/geno/fnf_geno_other/FNF_samples_with_replacements.txt"
    log:
        out="output/logs/categorize_femur_replicate_samples.out",
        err="output/logs/categorize_femur_replicate_samples.err"
    shell:
        """
        mkdir -p output/logs
        mkdir -p output/geno/pbs_geno_other
        mkdir -p output/geno/fnf_geno_other
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/categorize_samplesv3.sh 1> {log.out} 2> {log.err}
        
        # Move output files to expected locations
        mv /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/CTL_samples_with_replacements.txt {output.ctl_femur_replicate} 2>> {log.err}
        mv /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/FNF_samples_with_replacements.txt {output.fnf_femur_replicate} 2>> {log.err}
        """


rule reheader_subset:
    input:
        matched_pbs=rules.categorize_femur_replicate_samples.output.ctl_femur_replicate,
        matched_fnf=rules.categorize_femur_replicate_samples.output.fnf_femur_replicate
    output:
        pbs_vcf="output/geno/pbs_geno_other/pbs.rename_wCHR.vcf.gz",
        fnf_vcf="output/geno/fnf_geno_other/fnf.rename_wCHR.vcf.gz",
        pbs_temp="output/geno/pbs_geno_other/pbs.rename.vcf.gz",
        fnf_temp="output/geno/fnf_geno_other/fnf.rename.vcf.gz",
        mapping_pbs_txt='output/geno/pbs_geno_other/mappings/sample_name_mapping.txt',
        mapping_fnf_txt='output/geno/fnf_geno_other/mappings/sample_name_mapping.txt'
    params:
        vcf_copied=config['vcf_copy']
    log:
        pbsout="output/logs/PBS_rename_other.out",
        pbserr="output/logs/PBS_rename_other.err",
        fnfout="output/logs/FNF_rename_other.out",
        fnferr="output/logs/FNF_rename_other.err"
    shell:
        """
        
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/rename.sh {params.vcf_copied} {input.matched_pbs} {output.pbs_temp} {output.pbs_vcf} 1> {log.pbsout} 2> {log.pbserr}
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/rename.sh {params.vcf_copied} {input.matched_fnf} {output.fnf_temp} {output.fnf_vcf} 1> {log.fnfout} 2> {log.fnferr}
        """

rule reheader_imputed:
    input:
        imputed_vcf="/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/imputation/michigan_imputation_result/chr{chrom}.dose.vcf.gz",
        matched_pbs=rules.reheader_subset.output.mapping_pbs_txt,
        matched_fnf=rules.reheader_subset.output.mapping_fnf_txt
    output:
        pbs_vcf="output/geno/pbs_geno_other/wasp_vcf/chr{chrom}.pbs.rename.vcf.gz",
        fnf_vcf="output/geno/fnf_geno_other/wasp_vcf/chr{chrom}.fnf.rename.vcf.gz"
    params:
        chromosomes=list(range(1,23)),
        samtools_version = config['samtoolsVers']
    log:
        pbs="output/logs/imputed_PBS_OTHER_chr{chrom}_rename.log",
        fnf="output/logs/imputed_FNF_OTHER_chr{chrom}_rename.log"
    shell:
        """
        module load samtools/{params.samtools_version}
        # Get the chromosome number from wildcards
        chrom={wildcards.chrom}

        mkdir -p output/geno/pbs_geno_other/wasp_vcf output/geno/fnf_geno_other/wasp_vcf

        sed 's/^0_//' {input.matched_pbs} > output/geno/pbs_geno_other/mappings/imputed_pbs_{wildcards.chrom}_mapping.txt
        sed 's/^0_//' {input.matched_fnf} > output/geno/fnf_geno_other/mappings/imputed_fnf_{wildcards.chrom}_mapping.txt

        # PBS files
        bcftools reheader \
            --samples output/geno/pbs_geno_other/mappings/imputed_pbs_{wildcards.chrom}_mapping.txt \
            {input.imputed_vcf} \
            -o {output.pbs_vcf} 2> {log.pbs}

        # FNF files
        bcftools reheader \
            --samples output/geno/fnf_geno_other/mappings/imputed_fnf_{wildcards.chrom}_mapping.txt \
            {input.imputed_vcf} \
            -o {output.fnf_vcf} 2> {log.fnf}

        # Index files
        tabix -p vcf {output.pbs_vcf}
        tabix -p vcf {output.fnf_vcf}
        """

rule verifybamid:
    input:
        bam = rules.add_read_groups.output.bam,
        bai = rules.add_read_groups.output.bai,
        pbs_vcf = rules.reheader_subset.output.pbs_vcf,
        fnf_vcf = rules.reheader_subset.output.fnf_vcf
    output:
        selfSM = 'output/QC/{sampleName}_verifyBamID.selfSM',
        selfRG = 'output/QC/{sampleName}_verifyBamID.selfRG',
        bestRG = 'output/QC/{sampleName}_verifyBamID.bestRG',
        bestSM = 'output/QC/{sampleName}_verifyBamID.bestSM',
        depthRG = 'output/QC/{sampleName}_verifyBamID.depthRG',
        depthSM = 'output/QC/{sampleName}_verifyBamID.depthSM'
    params:
        verifybamid = config['verifybamid'],
        out_prefix = "output/QC/{sampleName}_verifyBamID"
    wildcard_constraints:
        sampleName = "|".join([
            "CQTL_AM7778FE_FNF_Femur_1",
            "CQTL_AM7754_FNF_Ankle_replicate", 
            "CQTL_AM7755_FNF_Ankle_replicate",
            "CQTL_AM7778FE_CTL_Femur_1",
            "CQTL_AM7754_CTL_Ankle_replicate",
            "CQTL_AM7755_CTL_Ankle_replicate"
        ])
    log:
        err = "output/logs/verifyBamID_{sampleName}.err"
    shell:
        """
        
        if [[ "{wildcards.sampleName}" == *"CTL"* ]]; then
            vcf_file={input.pbs_vcf}
        elif [[ "{wildcards.sampleName}" == *"FNF"* ]]; then
            vcf_file={input.fnf_vcf}
        else
            echo "Error: Sample type not recognized in {wildcards.sampleName}"
            exit 1
        fi
        
        echo "$vcf_file"

        {params.verifybamid} --vcf $vcf_file --bam {input.bam} --bai {input.bai} --best --out {params.out_prefix} 2> {log.err}
        """

rule extract_sample_ids:
    input:
        pbs_matched=rules.categorize_femur_replicate_samples.output.ctl_femur_replicate,
        fnf_matched=rules.categorize_femur_replicate_samples.output.fnf_femur_replicate
    output:
        sample_ids_pbs="output/wasp/pbs_subset_sample_ids.txt",
        sample_ids_fnf="output/wasp/fnf_subset_sample_ids.txt"
    log:
        pbs_log="output/logs/extract_sample_ids_pbs_subset.log",
        fnf_log="output/logs/extract_sample_ids_fnf_subset.log"
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
        pbs_haplotypes="output/wasp/pbs_rep/haplotypes.h5",
        pbs_snp_index="output/wasp/pbs_rep/snp_index.h5",
        pbs_snp_tab="output/wasp/pbs_rep/snp_tab.h5",
        fnf_haplotypes="output/wasp/fnf_rep/haplotypes.h5",
        fnf_snp_index="output/wasp/fnf_rep/snp_index.h5",
        fnf_snp_tab="output/wasp/fnf_rep/snp_tab.h5"
    params:
        wasp_path=config['wasp_dir'],
        chromsize=config['chrom_size'],
        wasp_ver=config['waspVer']
    log:
        pbs_err="output/logs/pbs_rep_snp2h5.err",
        fnf_err="output/logs/fnf_rep_snp2h5.err"
    shell:
        """
        module add wasp/{params.wasp_ver}
        mkdir -p output/wasp/pbs_rep output/wasp/fnf_rep

        # Process PBS VCF
        {params.wasp_path}/snp2h5/snp2h5 \
        --chrom {params.chromsize} \
        --format vcf \
        --haplotype {output.pbs_haplotypes} \
        --snp_index {output.pbs_snp_index} \
        --snp_tab {output.pbs_snp_tab} \
        output/geno/pbs_geno_other/wasp_vcf/chr*.pbs.rename.vcf.gz 2> {log.pbs_err}

        # Process FNF VCF
        {params.wasp_path}/snp2h5/snp2h5 \
        --chrom {params.chromsize} \
        --format vcf \
        --haplotype {output.fnf_haplotypes} \
        --snp_index {output.fnf_snp_index} \
        --snp_tab {output.fnf_snp_tab} \
        output/geno/fnf_geno_other/wasp_vcf/chr*.fnf.rename.vcf.gz 2> {log.fnf_err}
        """

rule find_intersecting_snps:
    input:
        bam=rules.add_read_groups.output.bam,
        pbs_snp_tab=rules.snp2h5.output.pbs_snp_tab,
        pbs_snp_index=rules.snp2h5.output.pbs_snp_index,
        pbs_haplotypes=rules.snp2h5.output.pbs_haplotypes,
        fnf_snp_tab=rules.snp2h5.output.fnf_snp_tab,
        fnf_snp_index=rules.snp2h5.output.fnf_snp_index,
        fnf_haplotypes=rules.snp2h5.output.fnf_haplotypes,
        sample_ids_pbs=rules.extract_sample_ids.output.sample_ids_pbs,
        sample_ids_fnf=rules.extract_sample_ids.output.sample_ids_fnf
    output:
        pbs_fq1="output/wasp/filter_remapped_reads/{sampleName}_RG.remap.fq1.gz",
        pbs_fq2="output/wasp/filter_remapped_reads/{sampleName}_RG.remap.fq2.gz",
        to_remap_bam="output/wasp/filter_remapped_reads/{sampleName}_RG.to.remap.bam",
        keep_bam="output/wasp/filter_remapped_reads/{sampleName}_RG.keep.bam"
    threads: 4
    wildcard_constraints:
        sampleName = "|".join([
            "CQTL_AM7778FE_FNF_Femur_1",
            "CQTL_AM7754_FNF_Ankle_replicate",
            "CQTL_AM7755_FNF_Ankle_replicate",
            "CQTL_AM7778FE_CTL_Femur_1",
            "CQTL_AM7754_CTL_Ankle_replicate",
            "CQTL_AM7755_CTL_Ankle_replicate"
        ])
    params:
        wasp_path=config['wasp_dir'],
        wasp_ver=config['waspVer']
    log:
        err="output/logs/find_intersecting_snps_pbs_rep_{sampleName}.err"
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
        rmdup_bam="output/wasp/rmdup/{sampleName}.rmdup.bam"
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

        """

rule rmblacklist_region:
    input:
        bam=rules.rmDupWASP.output.rmdup_bam
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
        
        samtools sort -o {output.sort_bam} {output.final_bam} 1> {log.out} 2> {log.err}
        samtools index {output.sort_bam} 1>> {log.out} 2>> {log.err}
        samtools flagstat {output.sort_bam} > {output.stats} 2>> {log.err}

        """


rule run_ataqv_subset:
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

rule macs3_frip_calc:
    input:
        saf="output/peaks/merged/allsamples_macs3_merged.saf",
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
                -a {input.saf} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}


        """
rule macs2_frip_calc:
    input:
        saf="output/peaks/merged/allsamples_macs2_merged.saf",
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
                -a {input.saf} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """
rule hmmratac_frip_calc:
    input:
        saf="output/peaks/merged/allsamples_hmmratac_merged.saf",
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

        #convert peak region to SAF
        #awk 'BEGIN{{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' indv.peakfile > indv.saf

        #count-indv
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {input.saf} \
                -o {output.count} \
                {input.bam} 2>> {log.err} 1>> {log.out}

        """

rule rocco_frip_calc:
    input:
        saf="output/peaks/merged/allsamples_rocco_merged.saf",
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

        #convert peak region to SAF
        #awk 'BEGIN{{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}} {{print "peak_"NR, $1, $2, $3, "."}}' input.peak > output.saf



        #count-indv
        featureCounts -p -F SAF \
                -T {threads} \
                --primary \
                -C \
                -a {input.saf} \
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


