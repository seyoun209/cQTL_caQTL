#!/usr/bin/env python3

import pandas as pd
import os, subprocess, glob, shutil

def OUT(*parts):
    return os.path.join(config["paths"]["genotype_outdir"], *[p for p in parts if p])

## Read in samplesheet
genos = pd.read_csv(config["geno"], sep=",")
genos = genos.astype(str)

## Concatenate Genotyping_Directory to File_Prefix for full path
genos['Name'] = genos[['Genotyping_Directory', 'File_Prefix']].apply(lambda row: os.path.join(*row), axis=1)

mergeName = '_'.join(genos.Proj.unique()) + '_' + '_'.join(genos['Batch'])

onsuccess:
    print("Genotype data processing finished successfully!")
    print(f"VCFs split by chromosome in '{OUT('imputation')}' are ready for imputation.")
    print(f"Ancestry plot can be found at {OUT('ancestry', mergeName + '_ancestry.pdf')}")
    
    ## Clean up intermediate files
    for f in glob.glob(f"{OUT('binary')}/*.hh"):
        os.remove(f)
    for f in glob.glob(f"{OUT('merge')}/*.hh"):
        os.remove(f)
    for f in glob.glob(f"{OUT('ancestry')}/*.hh"):
        os.remove(f)
    for f in glob.glob(f"{OUT('ancestry')}/*.nosex"):
        os.remove(f)

rule all:
    input:
        expand(OUT("imputation", f"{mergeName}_chr{{chr}}.recode.vcf.gz"), chr=range(1, 23)),
        OUT("reports", mergeName + ".sexcheck"),
        OUT("reports", mergeName + ".genome"),
        OUT("ancestry", mergeName + "_ref_merge.pca.evec"),
        OUT("ancestry", mergeName + "_ref_merge.eval"),
        OUT("ancestry", mergeName + "_ancestry.pdf")

rule convertBinary:
    input:
        expand("{name}{ext}", name=genos['Name'], ext=[".map",".ped"])
    output:
        temp(expand(OUT("binary", "{name}{ext}"), name=genos['Batch'], ext=[".bed",".bim",".fam"]))
    params:
        binary_dir=OUT("binary"),
        logs_dir=OUT("binary", "logs")
    run:
        os.makedirs(params.logs_dir, exist_ok=True)
        for row in genos.itertuples():
            # Use the exact same approach as your working version, just with module loading
            cmd = f"module load plink/{config['plink']} && plink --file {row.Name} --make-bed --allow-extra-chr --out {params.binary_dir}/{row.Batch} 1> {params.logs_dir}/{row.Batch}_binary.out"
            p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            shutil.move(f"{params.binary_dir}/{row.Batch}.log", f"{params.logs_dir}/{row.Batch}_binary.err")

rule removeSamples:
    input:
        rules.convertBinary.output
    output:
        temp(expand(OUT("binary", "{name}_filter{ext}"), name=genos['Batch'], ext=[".bed",".bim",".fam"])),
        removeFile=temp("remove.plink")
    params:
        removeSamples=config['remove'],
        binary_dir=OUT("binary"),
        logs_dir=OUT("binary", "logs")
    run:
        subprocess.call(['python3', 'scripts/geno_workflow/removeSamples.py', params.removeSamples])
        
        for row in genos.itertuples():
            cmd = f"module load plink/{config['plink']} && plink --bfile {params.binary_dir}/{row.Batch} --remove remove.plink --make-bed --out {params.binary_dir}/{row.Batch}_filter 1> {params.logs_dir}/{row.Batch}_filter.out"
            p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            shutil.move(f"{params.binary_dir}/{row.Batch}_filter.log", f"{params.logs_dir}/{row.Batch}_filter.err")

rule mergeData:
    input:
        rules.removeSamples.output
    output:
        temp("mergeList.txt"),
        bed=temp(OUT("merge", mergeName + ".bed")),
        bim=temp(OUT("merge", mergeName + ".bim")),
        fam=temp(OUT("merge", mergeName + ".fam"))
    params:
        file=genos.loc[0, "Batch"] + "_filter", 
        prefix=mergeName,
        binary_dir=OUT("binary"),
        merge_dir=OUT("merge"),
        logs_dir=OUT("merge", "logs"),
        reports_dir=OUT("reports")
    log:
        out=OUT("merge", "logs", mergeName + ".out")
    run:
        os.makedirs(params.logs_dir, exist_ok=True)
        os.makedirs(params.reports_dir, exist_ok=True)

        mergeList = open('mergeList.txt', 'w')
        for row in genos[1:].itertuples():
            mergeList.write(f"{params.binary_dir}/{row.Batch}_filter\n")
        mergeList.close()

        cmd = f"module load plink/{config['plink']} && plink --bfile {params.binary_dir}/{params.file} --merge-list mergeList.txt --merge-equal-pos --make-bed --out {params.merge_dir}/{params.prefix} 1> {log.out}"
        p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        i = 1
        while p.returncode != 0:
            grep = f"cat {params.merge_dir}/{mergeName}.log | tr '\\r\\n' ' ' | grep -o 'Error: --merge-equal-pos failure.  Variants .*' | awk '{{print $5}}' | tr -d \"'\" >> mergeFails.txt && cat {params.merge_dir}/{mergeName}.log | tr '\\r\\n' ' ' | grep -o 'Error: --merge-equal-pos failure.  Variants .*' | awk '{{print $7}}' | tr -d \"'\" >> mergeFails.txt"
            p1 = subprocess.run(grep, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            mergeList = open('mergeList.txt', 'w')
            for row in genos.itertuples():
                if i != 1:
                    for f in glob.glob(f"{params.binary_dir}/{row.Batch}_filter_{i-1}.*"):
                        os.remove(f)
                plink_subset = f"plink --bfile {params.binary_dir}/{row.Batch}_filter --exclude mergeFails.txt --make-bed --out {params.binary_dir}/{row.Batch}_filter_{i} 1> {params.logs_dir}/resubset.out"
                p2 = subprocess.run(plink_subset, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if row.Index > 0:
                    mergeList.write(f"{params.binary_dir}/{row.Batch}_filter_{i}\n")
            mergeList.close()

            remerge = f"plink --bfile {params.binary_dir}/{params.file}_{i} --merge-list mergeList.txt --merge-equal-pos --make-bed --out {params.merge_dir}/{params.prefix} 1>{log.out}"   
            p = subprocess.run(remerge, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            i += 1
            
        shutil.move("mergeFails.txt", f"{params.reports_dir}/{mergeName}_mergeFails.txt")
        shutil.move(f"{params.merge_dir}/{mergeName}.log", f"{params.logs_dir}/{mergeName}.err")
        if os.path.exists(f"{params.logs_dir}/resubset.out"):
            os.remove(f"{params.logs_dir}/resubset.out")

rule qc:
    input:
        rules.mergeData.output
    output:
        bed=temp(OUT("merge", mergeName + "_qc.bed")),
        bim=temp(OUT("merge", mergeName + "_qc.bim")),
        fam=temp(OUT("merge", mergeName + "_qc.fam"))
    params:
        prefix=OUT("merge", mergeName),
        outpre=mergeName,
        geno=config['miss_call'],
        hwe=config['hwe'],
        maf=config['maf'],
        merge_dir=OUT("merge"),
        logs_dir=OUT("merge", "logs")
    log:
        out=OUT("merge", "logs", mergeName + "_qc.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --geno {params.geno} --hwe {params.hwe} --maf {params.maf} --make-bed --out {params.merge_dir}/{params.outpre}_qc 1> {log.out}
        mv {params.merge_dir}/{params.outpre}_qc.log {params.logs_dir}/{params.outpre}_qc.err
        """

rule sexCheck:
    input: 
        rules.qc.output
    output:
        OUT("reports", mergeName + ".sexcheck")
    params:
        prefix=OUT("merge", mergeName + "_qc"),
        outpre=mergeName,
        reports_dir=OUT("reports")
    shell:
        """
        mkdir -p {params.reports_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --check-sex --out {params.reports_dir}/{params.outpre}

        if [[ -f {params.reports_dir}/{params.outpre}.hh ]]; then
            rm {params.reports_dir}/{params.outpre}.hh
        fi
        rm {params.reports_dir}/{params.outpre}.log
        """

rule autosomePrune:
    input:
        rules.qc.output
    output:
        bed=OUT("merge", mergeName + "_qcautosome.bed"),
        bim=OUT("merge", mergeName + "_qcautosome.bim"),
        fam=OUT("merge", mergeName + "_qcautosome.fam")
    params:
        prefix=OUT("merge", mergeName + "_qc"),
        outpre=mergeName,
        merge_dir=OUT("merge"),
        logs_dir=OUT("merge", "logs")
    log:
        out=OUT("merge", "logs", mergeName + "_qcautosome.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --autosome --make-bed --out {params.merge_dir}/{params.outpre}_qcautosome 1> {log.out}
        mv {params.merge_dir}/{params.outpre}_qcautosome.log {params.logs_dir}/{params.outpre}_qcautosome.err
        """

rule snpFlip:
    input:
        rules.autosomePrune.output
    output:
        rev=temp(OUT("imputation", "snpflip.reverse")),
        amb=temp(OUT("imputation", "snpflip.ambiguous")),
        an_bim=temp(OUT("imputation", "snpflip.annotated_bim"))
    params:
        prefix=OUT("merge", mergeName + "_qcautosome"),
        sequence=config['sequence'],
        imputation_dir=OUT("imputation"),
        python_ver=config['python']
    shell:
        """
        module load python/{params.python_ver}
        mkdir -p {params.imputation_dir}
        
        
        snpflip --fasta-genome={params.sequence} --bim={params.prefix}.bim --output-prefix={params.imputation_dir}/snpflip
        """

rule flipSnps:
    input:
        rules.autosomePrune.output, 
        rules.snpFlip.output.rev
    output:
        bed=temp(OUT("imputation", mergeName + "_qcautosome_snpflip.bed")),
        bim=temp(OUT("imputation", mergeName + "_qcautosome_snpflip.bim")),
        fam=temp(OUT("imputation", mergeName + "_qcautosome_snpflip.fam"))
    params:
        prefix=OUT("merge", mergeName + "_qcautosome"),
        outpre=mergeName,
        imputation_dir=OUT("imputation"),
        logs_dir=OUT("imputation", "logs")
    log:
        out=OUT("imputation", "logs", mergeName + "_qcautosome_snpflip.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --flip {input[1]} --make-bed --out {params.imputation_dir}/{params.outpre}_qcautosome_snpflip 1> {log.out}
        mv {params.imputation_dir}/{params.outpre}_qcautosome_snpflip.log {params.logs_dir}/{params.outpre}_qcautosome_snpflip.err 
        """

rule convertVCF:
    input:
        rules.flipSnps.output
    output:
        temp(OUT("imputation", mergeName + ".vcf"))
    params:
        prefix=OUT("imputation", mergeName + "_qcautosome_snpflip"),
        outpre=mergeName,
        imputation_dir=OUT("imputation"),
        logs_dir=OUT("imputation", "logs")
    log:
        out=OUT("imputation", "logs", mergeName + "_vcf.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --recode vcf --out {params.imputation_dir}/{params.outpre} 1> {log.out}
        mv {params.imputation_dir}/{params.outpre}.log {params.logs_dir}/{params.outpre}_vcf.err
        """

rule splitVCF:
    input:
        rules.convertVCF.output
    output:
        temp(OUT("imputation", mergeName + ".vcf.gz")),
        temp(OUT("imputation", mergeName + ".vcf.gz.tbi")),
        expand(OUT("imputation", f"{mergeName}_chr{{chr}}.recode.vcf.gz"), chr=range(1, 23))
    params:
        prefix=OUT("imputation", mergeName),
        samtools_ver=config['samtools']
    shell:
        """
        module load samtools/{params.samtools_ver}
        module load vcftools

        bgzip {input}
        tabix -p vcf {params.prefix}.vcf.gz
        
        for chr in {{1..22}}; do
            vcftools --gzvcf {params.prefix}.vcf.gz --chr ${{chr}} --recode --recode-INFO-all --out {params.prefix}_chr${{chr}}
            bgzip {params.prefix}_chr${{chr}}.recode.vcf
        done
        """

rule ldPrune:
    input:
        rules.qc.output
    output:
        infile=temp(OUT("ancestry", mergeName + "_ldprune.prune.in")),
        outfile=temp(OUT("ancestry", mergeName + "_ldprune.prune.out"))
    params:
        prefix=OUT("merge", mergeName + "_qc"),
        outpre=mergeName,
        ancestry_dir=OUT("ancestry"),
        logs_dir=OUT("ancestry", "logs")
    log:
        out=OUT("ancestry", "logs", mergeName + "_ldprune.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --indep-pairwise 50 5 0.1 --out {params.ancestry_dir}/{params.outpre}_ldprune 1> {log.out}
        mv {params.ancestry_dir}/{params.outpre}_ldprune.log {params.logs_dir}/{params.outpre}_ldprune.err
        """

rule pruneData:
    input:
        rules.qc.output,
        ldList=rules.ldPrune.output.infile
    output:
        bed=temp(OUT("ancestry", mergeName + "_ldautosome.bed")),
        bim=temp(OUT("ancestry", mergeName + "_ldautosome.bim")),
        fam=temp(OUT("ancestry", mergeName + "_ldautosome.fam"))
    params:
        prefix=OUT("merge", mergeName + "_qc"),
        outpre=mergeName,
        ancestry_dir=OUT("ancestry"),
        logs_dir=OUT("ancestry", "logs")
    log:
        out=OUT("ancestry", "logs", mergeName + "_ldautosome.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --extract {input.ldList} --autosome --make-bed --out {params.ancestry_dir}/{params.outpre}_ldautosome 1> {log.out}
        mv {params.ancestry_dir}/{params.outpre}_ldautosome.log {params.logs_dir}/{params.outpre}_ldautosome.err
        """

rule pairwiseIBD:
    input: 
        rules.pruneData.output
    output:
        OUT("reports", mergeName + ".genome")
    params:
        prefix=OUT("ancestry", mergeName + "_ldautosome"),
        outpre=mergeName,
        reports_dir=OUT("reports")
    shell:
        """
        mkdir -p {params.reports_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.prefix} --genome --out {params.reports_dir}/{params.outpre}
        rm {params.reports_dir}/{params.outpre}.log
        """

rule pruneReference:
    input:
        rules.pruneData.output
    output:
        bed=temp(OUT("ancestry", "refPrune.bed")),
        bim=temp(OUT("ancestry", "refPrune.bim")),
        fam=temp(OUT("ancestry", "refPrune.fam"))
    params:
        prefix=OUT("ancestry", mergeName + "_ldautosome"),
        outpre=mergeName,
        ref=config['ref'],
        ancestry_dir=OUT("ancestry"),
        logs_dir=OUT("ancestry", "logs")
    log:
        out=OUT("ancestry", "logs", "refPrune.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}

        awk '{{print $1" "$4" "$4" "$2}}' {params.prefix}.bim > {params.ancestry_dir}/{params.outpre}_list.txt

        plink --bfile {params.ref} --extract range {params.ancestry_dir}/{params.outpre}_list.txt --make-bed --out {params.ancestry_dir}/refPrune 1> {log.out}
        mv {params.ancestry_dir}/refPrune.log {params.logs_dir}/refPrune.err

        rm {params.ancestry_dir}/{params.outpre}_list.txt
        """

rule setA1:
    input:
        rules.pruneData.output,
        rules.pruneReference.output
    output:
        bed=temp(OUT("ancestry", mergeName + "_ldautosome_A1.bed")),
        bim=temp(OUT("ancestry", mergeName + "_ldautosome_A1.bim")),
        fam=temp(OUT("ancestry", mergeName + "_ldautosome_A1.fam"))
    params:
        outpre=mergeName,
        ancestry_dir=OUT("ancestry"),
        logs_dir=OUT("ancestry", "logs"),
        python_ver=config['python']
    log:
        out=OUT("ancestry", "logs", mergeName + "_ldautosome_A1.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load python/{params.python_ver}
        module load plink/{config[plink]}
        
        awk '{{print $1":"$4" "$5}}' {params.ancestry_dir}/refPrune.bim > {params.ancestry_dir}/setA1.txt
        
        python3 scripts/geno_workflow/join.py {params.ancestry_dir}/{params.outpre}_ldautosome.bim {params.ancestry_dir}/setA1.txt {params.ancestry_dir}/A1.txt
        
        plink --bfile {params.ancestry_dir}/{params.outpre}_ldautosome --a1-allele {params.ancestry_dir}/A1.txt --make-bed --out {params.ancestry_dir}/{params.outpre}_ldautosome_A1 1> {log.out}
        
        mv {params.ancestry_dir}/{params.outpre}_ldautosome_A1.log {params.logs_dir}/{params.outpre}_ldautosome_A1.err

        rm {params.ancestry_dir}/setA1.txt {params.ancestry_dir}/A1.txt
        """

rule mismatchAlleles:
    input:
        rules.pruneReference.output,
        rules.setA1.output.bed, rules.setA1.output.bim, rules.setA1.output.fam
    output:
       temp(expand(OUT("ancestry", "{name}{ext}"), name=[mergeName + "_ldautosome_A1_final", "refPrune_final"], ext=[".bim", ".bed", ".fam"]))
    params:
        outpre=mergeName,
        ancestry_dir=OUT("ancestry"),
        logs_dir=OUT("ancestry", "logs"),
        r_ver=config.get('r', '4.5.0')
    log:
        out1=OUT("ancestry", "logs", mergeName + "_mismatch.out"),
        out2=OUT("ancestry", "logs", "refPrune_mismatch.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load r/{params.r_ver}
        module load plink/{config[plink]}
        
        Rscript scripts/geno_workflow/mismatchAlleles.r {params.ancestry_dir}/{params.outpre}_ldautosome_A1.bim {params.ancestry_dir}/refPrune.bim

        plink --bfile {params.ancestry_dir}/{params.outpre}_ldautosome_A1 --exclude {params.ancestry_dir}/data_missnp.txt --make-bed --out {params.ancestry_dir}/{params.outpre}_ldautosome_A1_final 1> {log.out1}
        plink --bfile {params.ancestry_dir}/refPrune --exclude {params.ancestry_dir}/ref_missnp.txt --make-bed --out {params.ancestry_dir}/refPrune_final 1> {log.out2}

        mv {params.ancestry_dir}/{params.outpre}_ldautosome_A1_final.log {params.logs_dir}/{params.outpre}_ldautosome_A1_final.err
        mv {params.ancestry_dir}/refPrune_final.log {params.logs_dir}/refPrune_final.err 

        rm {params.ancestry_dir}/data_missnp.txt {params.ancestry_dir}/ref_missnp.txt
        """

rule mergeDataReference:
    input:
        rules.mismatchAlleles.output
    output:
        ped=OUT("ancestry", mergeName + "_ref_merge.ped"),
        mapf=OUT("ancestry", mergeName + "_ref_merge.map"),
    params:
        prefix=OUT("ancestry", mergeName + "_ldautosome_A1_final"),
        outpre=mergeName,
        ancestry_dir=OUT("ancestry"),
        logs_dir=OUT("ancestry", "logs")
    log:
        out=OUT("ancestry", "logs", mergeName + "_ref_merge.out")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load plink/{config[plink]}
        
        plink --bfile {params.ancestry_dir}/refPrune_final --bmerge {params.prefix}.bed {params.prefix}.bim {params.prefix}.fam --merge-equal-pos --recode --out {params.ancestry_dir}/{params.outpre}_ref_merge 1> {log.out}
        mv {params.ancestry_dir}/{params.outpre}_ref_merge.log {params.logs_dir}/{params.outpre}_ref_merge.err

        awk '{{print $1, substr($2,1,39), "0", $4}}' {params.ancestry_dir}/{params.outpre}_ref_merge.map > {params.ancestry_dir}/{params.outpre}_ref_merge_temp.map
        rm {params.ancestry_dir}/{params.outpre}_ref_merge.map
        mv {params.ancestry_dir}/{params.outpre}_ref_merge_temp.map {params.ancestry_dir}/{params.outpre}_ref_merge.map

        awk '{{$6=1 ; print ;}}' {params.ancestry_dir}/{params.outpre}_ref_merge.ped > {params.ancestry_dir}/{params.outpre}_ref_merge_temp.ped
        rm {params.ancestry_dir}/{params.outpre}_ref_merge.ped
        mv {params.ancestry_dir}/{params.outpre}_ref_merge_temp.ped {params.ancestry_dir}/{params.outpre}_ref_merge.ped
        """
        
rule convertEigenstrat:
    input:
        ped=rules.mergeDataReference.output.ped,
        mapf=rules.mergeDataReference.output.mapf
    output:
        OUT("ancestry", mergeName + "_ref_merge.eigenstratgeno"),
        OUT("ancestry", mergeName + "_ref_merge.snp"),
        OUT("ancestry", mergeName + "_ref_merge.ind")
    params:
        prefix=OUT("ancestry", mergeName + "_ref_merge"),
        eigensoft=config['eigensoft'],
        ancestry_dir=OUT("ancestry")
    shell:
        """
        echo "genotypename: "{input.ped} >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "snpname: "{input.mapf} >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "indivname: "{input.ped} >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "outputformat: EIGENSTRAT" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "genotypeoutname: "{params.prefix}".eigenstratgeno" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "snpoutname: "{params.prefix}".snp" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "indivoutname: "{params.prefix}".ind" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "familynames: NO" >> {params.ancestry_dir}/par.PED.EIGENSTRAT

        {params.eigensoft}/convertf -p {params.ancestry_dir}/par.PED.EIGENSTRAT

        rm {params.ancestry_dir}/par.PED.EIGENSTRAT
        """

rule runEigenstrat:
    input:
        rules.convertEigenstrat.output
    output:
        OUT("ancestry", mergeName + "_ref_merge.pca.evec"),
        OUT("ancestry", mergeName + "_ref_merge.eval")
    params:
        prefix=OUT("ancestry", mergeName + "_ref_merge"),
        outpre=mergeName,
        eigensoft=config['eigensoft'],
        logs_dir=OUT("ancestry", "logs")
    shell:
        """
        mkdir -p {params.logs_dir}
        
        {params.eigensoft}/smartpca.perl -i {params.prefix}.eigenstratgeno -a {params.prefix}.snp -b {params.prefix}.ind -o {params.prefix}.pca -p {params.prefix}.plot -e {params.prefix}.eval -l {params.prefix}.log -m 0

        mv {params.prefix}.log {params.logs_dir}/{params.outpre}_smartpca.err

        sed '1d' {params.prefix}.pca.evec > tempfile
        rm {params.prefix}.pca.evec
        mv tempfile {params.prefix}.pca.evec
        """

rule plotAncestry:
    input:
        rules.runEigenstrat.output
    output:
        plot=OUT("ancestry", mergeName + "_ancestry.pdf"),
        predictions=OUT("ancestry", mergeName + "_predictedAncestry.csv")
    params:
        panel=config['panel'],
        pca=OUT("ancestry", mergeName + "_ref_merge.pca.evec"),
        name=config['pop_name'],
        r_ver=config.get('r', '4.5.0'),
        ancestry_dir=OUT("ancestry")
    shell:
        """
        mkdir -p {params.ancestry_dir}
        module load r/{params.r_ver}
        Rscript scripts/geno_workflow/getAncestry.r {params.pca} {params.panel} {params.name} {output.plot} {output.predictions}
        """
