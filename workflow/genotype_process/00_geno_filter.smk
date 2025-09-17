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
        expand(OUT("imputation", f"{mergeName}_chr{{chr}}.vcf.gz"), chr=range(1, 23)),
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
        logs_dir=OUT("binary", "logs"),
        plink_ver=config['plink']
    run:
        os.makedirs(params.logs_dir, exist_ok=True)
        for row in genos.itertuples():
            cmd = f"module load plink/{params.plink_ver} && plink --file {row.Name} --make-bed --allow-extra-chr --out {params.binary_dir}/{row.Batch}"
            p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if os.path.exists(f"{params.binary_dir}/{row.Batch}.log"):
                shutil.move(f"{params.binary_dir}/{row.Batch}.log", f"{params.logs_dir}/{row.Batch}_binary.err")

rule removeSamples:
    input:
        rules.convertBinary.output
    output:
        temp(expand(OUT("binary", "{name}_clean{ext}"), name=genos['Batch'], ext=[".bed",".bim",".fam"])),
        removeFile=temp("remove.plink")
    params:
        removeSamples=config['remove'],
        binary_dir=OUT("binary"),
        logs_dir=OUT("binary", "logs"),
        plink_ver=config['plink']
    run:
        import os
        import subprocess
        import shutil
        
        os.makedirs(params.logs_dir, exist_ok=True)
        
        # Generate remove.plink file
        result = subprocess.run(['python3', 'scripts/geno_workflow/removeSamples.py', params.removeSamples], 
                              capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"removeSamples.py failed: {result.stderr}")
            # Create empty remove.plink on failure
            with open("remove.plink", "w") as f:
                pass
        
        # Check if there are samples to remove
        has_samples_to_remove = os.path.exists("remove.plink") and bool(open("remove.plink").read().strip())

        # Process each batch
        for row in genos.itertuples():
            batch_base = f"{params.binary_dir}/{row.Batch}"
            clean_file = f"{batch_base}_clean"

            # Remove samples and clean variant IDs in one step
            if has_samples_to_remove:
                cmd = f"module load plink/{params.plink_ver} && plink --bfile {batch_base} --remove remove.plink --set-missing-var-ids @:# --make-bed --out {clean_file}"
            else:
                cmd = f"module load plink/{params.plink_ver} && plink --bfile {batch_base} --set-missing-var-ids @:# --make-bed --out {clean_file}"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Error processing {row.Batch}: {result.stderr}")
                raise RuntimeError(f"PLINK processing failed for {row.Batch}")
            
            # Move log file
            if os.path.exists(f"{clean_file}.log"):
                shutil.move(f"{clean_file}.log", f"{params.logs_dir}/{row.Batch}_clean.log")

rule mergeData:
    input:
        rules.removeSamples.output
    output:
        bed=OUT("merge", mergeName + ".bed"),
        bim=OUT("merge", mergeName + ".bim"),
        fam=OUT("merge", mergeName + ".fam")
    params:
        file=genos.loc[0, "Batch"] + "_clean",
        prefix=mergeName,
        merge_dir=OUT("merge"),
        binary_dir=OUT("binary"),
        reports_dir=OUT("reports"),
        plink_ver=config['plink']
    log:
        out=OUT("merge", "logs", mergeName + ".out")
    run:
        os.makedirs(f"{params.merge_dir}/logs", exist_ok=True)
        os.makedirs(params.reports_dir, exist_ok=True)

        # If only one batch, just copy files
        if len(genos) == 1:
            batch = genos.loc[0, "Batch"]
            for ext in [".bed", ".bim", ".fam"]:
                shutil.copy2(f"{params.binary_dir}/{batch}_clean{ext}",
                           f"{params.merge_dir}/{params.prefix}{ext}")
        else:
            # Multiple batches - simple merge
            mergeList = open('mergeList.txt', 'w')
            for row in genos[1:].itertuples():
                mergeList.write(f"{params.binary_dir}/" + row.Batch + "_clean\n")
            mergeList.close()

            # Try merge - let PLINK handle conflicts automatically
            cmd = f"module load plink/{params.plink_ver} && plink --bfile {params.binary_dir}/" + params.file + f" --merge-list mergeList.txt --make-bed --out {params.merge_dir}/" + params.prefix
            p = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if p.returncode != 0:
                print("Merge failed:", p.stderr)
                raise RuntimeError("PLINK merge failed")

            os.remove("mergeList.txt")

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

rule create_frequency_file:
    input:
        rules.autosomePrune.output
    output:
        freq=OUT("imputation", mergeName + "_qcautosome.frq")
    params:
        prefix=OUT("merge", mergeName + "_qcautosome"),
        outpre=mergeName,
        imputation_dir=OUT("imputation"),
        plink_ver=config['plink']
    shell:
        """
        mkdir -p {params.imputation_dir}
        module load plink/{params.plink_ver}
        
        plink --bfile {params.prefix} --freq --out {params.imputation_dir}/{params.outpre}_qcautosome
        """

rule allele_swap_check:
    input:
        bed=rules.autosomePrune.output.bed,
        bim=rules.autosomePrune.output.bim,
        fam=rules.autosomePrune.output.fam,
        freq=rules.create_frequency_file.output.freq
    output:
        bed=OUT("imputation", "allele_swap.bed"),
        bim=OUT("imputation", "allele_swap.bim"),
        fam=OUT("imputation", "allele_swap.fam")
    params:
        input_prefix=OUT("merge", mergeName + "_qcautosome"),
        merge_dir=OUT("merge"),
        plink_ver=config['plink'],
        perl_ver=config['perl'],
        rayner_script=config['rayner_script'],
        rayner_legend=config['rayner_legend'],
        ref_path=config['ref']
    shell:
        """
        module load plink/{params.plink_ver}
        module load perl/{params.perl_ver}
        
        # Run right  reference genome
        if [[ "{params.ref_path}" == *"GRCh38"* ]]; then
            echo "Detected GRCh38 - running Rayner with -h flag"
            perl {params.rayner_script} -b {params.input_prefix}.bim -f {input.freq} -r {params.rayner_legend} -h
        else
            echo "Detected GRCh37/hg19 - running Rayner with -g -p ALL flags"
            perl {params.rayner_script} -b {params.input_prefix}.bim -f {input.freq} -r {params.rayner_legend} -g -p ALL
        fi 

        # Execute fixes
        cd {params.merge_dir} && bash Run-plink.sh && cd -
        
        # Copy corrected PLINK files
        cp {params.merge_dir}/CQTL_COA9_10_qcautosome-updated.bed {output.bed}
        cp {params.merge_dir}/CQTL_COA9_10_qcautosome-updated.bim {output.bim}
        cp {params.merge_dir}/CQTL_COA9_10_qcautosome-updated.fam {output.fam}
       
        """

rule vcf_zip:
    input:
        bed=rules.allele_swap_check.output.bed
    output:
        expand(OUT("imputation", f"{mergeName}_chr{{chr}}.vcf.gz"), chr=range(1, 23)),
    params:
        merge_dir=OUT("merge"),
        output_prefix=OUT("imputation", mergeName),
        samtools_ver=config['samtools'],
        ref_path=config['ref'],
        add_chr_file=config['add_chr']
    shell:
        """
        module load samtools/{params.samtools_ver}
        
        for chr in {{1..22}}; do
            source_file="{params.merge_dir}/CQTL_COA9_10_qcautosome-updated-chr${{chr}}.vcf"
            target_file="{params.output_prefix}_chr${{chr}}.vcf"
            
            # Add chr prefix if hg38, otherwise use as-is
            if [[ "{params.ref_path}" == *"GRCh38"* ]]; then
                bcftools annotate --rename-chrs {params.add_chr_file} "$source_file" -o "$target_file"
            else
                mv "$source_file" "$target_file"
            fi
            
            bgzip "$target_file"
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
