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
        expand(OUT("imputation", f"{mergeName}_chr{{chr}}_hg38.vcf.gz"), chr=range(1, 23)),
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
        logs_dir=OUT("binary", "logs"),
        plink_ver=config['plink']
    run:
        # Generate remove.plink file (empty if no samples to remove)
        subprocess.call(['python3', 'scripts/geno_workflow/removeSamples.py', params.removeSamples])

        # Check if there are samples to remove
        has_samples_to_remove = os.path.exists("remove.plink") and bool(open("remove.plink").read().strip())

        # Process each batch
        for row in genos.itertuples():
            batch_base = f"{params.binary_dir}/{row.Batch}"
            filter_file = f"{batch_base}_filter"

            # Step 1: Remove samples if needed
            if has_samples_to_remove:
                cmd = f"module load plink/{params.plink_ver} && plink --bfile {batch_base} --remove remove.plink --make-bed --out {filter_file}"
            else:
                cmd = f"module load plink/{params.plink_ver} && plink --bfile {batch_base} --make-bed --out {filter_file}"
            subprocess.run(cmd, shell=True, check=True)

            # Step 2: Clean variant IDs for merging
            clean_file = f"{params.binary_dir}/{row.Batch}_clean"
            clean_cmd = f"module load plink/{params.plink_ver} && plink --bfile {filter_file} --set-missing-var-ids @:# --make-bed --out {clean_file}"
            subprocess.run(clean_cmd, shell=True, check=True)

rule mergeData:
    input:
        geno_csv="geno.csv"
    output:
        bed=f"{params.merge_dir}/{params.prefix}.bed",
        bim=f"{params.merge_dir}/{params.prefix}.bim",
        fam=f"{params.merge_dir}/{params.prefix}.fam"
    params:
        merge_dir=params.merge_dir,
        binary_dir=params.binary_dir,
        plink_ver=params.plink_ver,
        prefix=params.prefix
    run:
        import pandas as pd, os, subprocess, shutil

        # Load batches
        genos = pd.read_csv(input.geno_csv)

        os.makedirs(params.merge_dir, exist_ok=True)

        # Only one batch → just copy
        if len(genos) == 1:
            batch = genos.loc[0, "Batch"]
            for ext in [".bed", ".bim", ".fam"]:
                shutil.copy2(f"{params.binary_dir}/{batch}_clean{ext}",
                             f"{params.merge_dir}/{params.prefix}{ext}")
            with open("mergeList.txt", "w") as f:
                pass
        else:
            # Write merge list
            with open("mergeList.txt", "w") as f:
                for row in genos[1:].itertuples():
                    f.write(f"{params.binary_dir}/{row.Batch}_clean\n")

            # Run PLINK merge
            first_batch = f"{params.binary_dir}/{genos.loc[0, 'Batch']}_clean"
            merged_out = f"{params.merge_dir}/{params.prefix}"

            cmd = (
                f"module load plink/{params.plink_ver} && "
                f"plink --bfile {first_batch} "
                f"--merge-list mergeList.txt "
                f"--merge-equal-pos --make-bed "
                f"--out {merged_out}"
            )
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                print(result.stderr)
                raise RuntimeError("PLINK merge failed")

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
        bed=rules.autosomePrune.output.bed,
        bim=rules.autosomePrune.output.bim,
        fam=rules.autosomePrune.output.fam,
        flip_file=rules.snpFlip.output.rev
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
        
        plink --bfile {params.prefix} --flip {input.flip_file} --make-bed --out {params.imputation_dir}/{params.outpre}_qcautosome_snpflip 1> {log.out}
        
        if [[ -f {params.imputation_dir}/{params.outpre}_qcautosome_snpflip.log ]]; then
            mv {params.imputation_dir}/{params.outpre}_qcautosome_snpflip.log {params.logs_dir}/{params.outpre}_qcautosome_snpflip.err
        fi
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

rule addChrPrefix:
    input:
        vcf=rules.convertVCF.output
    output:
        temp(OUT("imputation", mergeName + "_chr.vcf.gz"))
    params:
        temp_vcf=OUT("imputation", mergeName + "_chr.vcf"),
        mapping_file=config['add_chr'],
        samtools_ver=config['samtools']
    log:
        out=OUT("imputation", "logs", mergeName + "_addchr.out")
    shell:
        """
        module load samtools/{params.samtools_ver}

        # Use bcftools to rename chromosomes with your existing mapping file
        bcftools annotate --rename-chrs {params.mapping_file} {input.vcf} -o {params.temp_vcf} 2> {log.out}

        # Compress the output
        bgzip {params.temp_vcf}

        echo "Successfully added chr prefix using bcftools with addchr.txt" >> {log.out}
        """

rule splitVCF:
    input:
        rules.addChrPrefix.output
    output:
        temp(OUT("imputation", mergeName + "_chr.vcf.gz.tbi")),
        expand(OUT("imputation", f"{mergeName}_chr{{chr}}_hg38.vcf.gz"), chr=range(1, 23))
    params:
        input_vcf=OUT("imputation", mergeName + "_chr.vcf.gz"),
        output_prefix=OUT("imputation", mergeName),
        samtools_ver=config['samtools']
    shell:
        """
        module load samtools/{params.samtools_ver}
        module load vcftools

        tabix -p vcf {params.input_vcf}
        
        for chr in {{1..22}}; do
            vcftools --gzvcf {params.input_vcf} --chr chr${{chr}} --recode --recode-INFO-all --out {params.output_prefix}_chr${{chr}}_hg38_temp
            mv {params.output_prefix}_chr${{chr}}_hg38_temp.recode.vcf {params.output_prefix}_chr${{chr}}_hg38.vcf
            bgzip {params.output_prefix}_chr${{chr}}_hg38.vcf
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
        r_ver=config['r']
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
        eig_ver=config['eigensoft'],
        ancestry_dir=OUT("ancestry")
    shell:
        """
        module load eigensoft/{params.eig_ver}
        
        echo "genotypename: {input.ped}" > {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "snpname: {input.mapf}" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "indivname: {input.ped}" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "outputformat: EIGENSTRAT" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "genotypeoutname: {params.prefix}.eigenstratgeno" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "snpoutname: {params.prefix}.snp" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "indivoutname: {params.prefix}.ind" >> {params.ancestry_dir}/par.PED.EIGENSTRAT
        echo "familynames: NO" >> {params.ancestry_dir}/par.PED.EIGENSTRAT

        convertf -p {params.ancestry_dir}/par.PED.EIGENSTRAT

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
        eig_ver=config['eigensoft'],
        logs_dir=OUT("ancestry", "logs")
    shell:
        """
        mkdir -p {params.logs_dir}
        module load eigensoft/{params.eig_ver}
        

        smartpca.perl -i {params.prefix}.eigenstratgeno -a {params.prefix}.snp -b {params.prefix}.ind -o {params.prefix}.pca -p {params.prefix}.plot -e {params.prefix}.eval -l {params.prefix}.log -m 0

        if [[ -f {params.prefix}.log ]]; then
            mv {params.prefix}.log {params.logs_dir}/{params.outpre}_smartpca.err
        fi

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
        r_ver=config['r'],
        ancestry_dir=OUT("ancestry")
    shell:
        """
        mkdir -p {params.ancestry_dir}
        module load r/{params.r_ver}
        Rscript scripts/geno_workflow/getAncestry.r {params.pca} {params.panel} {params.name} {output.plot} {output.predictions}
        """
