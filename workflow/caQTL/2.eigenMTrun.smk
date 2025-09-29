#!/usr/bin/env python3
import os

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

# Load analysis parameters from config
CONDITIONS = config['conditions']
VCF_GROUPS = config['vcf_groups']
CHROMOSOMES = config['chromosomes']
PC_VALUES = config['atac_pcs']
WINDOW_KB = config['window']

# Create condition to VCF group mapping
CONDITION_TO_GROUP = dict(zip(CONDITIONS, VCF_GROUPS))

onsuccess:
    print("EigenMT input preparation completed successfully!")
    print(f"Inputs available in: {OUT('caQTL', 'eigenMT', 'inputs')}")

rule all:
    input:
        # EigenMT inputs - 4 files per condition/PC/chromosome combination
        expand(OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "gen.positions_chr{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        expand(OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "genotypes_chr{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        expand(OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "phe.positions_chr{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        expand(OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "qtl_chr{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        #eigenMT results
        expand(OUT("caQTL", "eigenMT", "results", "window_{window}kb","{condition}", "pc{pc}",
                   "eigenMT_chr{chr}_{condition}_pc{pc}.txt"),condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB)

rule prepare_eigenmt_inputs:
    input:
        adjusted_results=OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                           "{condition}_pc{pc}_{window}kb_adjusted.txt"),
        traw=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "04_recode",
                           "{chr}." + CONDITION_TO_GROUP[wc.condition] + ".traw")
    output:
        gen_positions=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb", 
                         "{condition}", "pc{pc}",
                         "gen.positions_chr{chr}_{condition}_pc{pc}.txt"),
        genotypes=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                     "{condition}", "pc{pc}",
                     "genotypes_chr{chr}_{condition}_pc{pc}.txt"),
        phe_positions=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                         "{condition}", "pc{pc}",
                         "phe.positions_chr{chr}_{condition}_pc{pc}.txt"),
        qtl=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
               "{condition}", "pc{pc}",
               "qtl_chr{chr}_{condition}_pc{pc}.txt")
    params:
        r_version=config['r']
    threads: 2
    resources:
        mem_mb=16000,
        time="1:00:00"
    log:
        OUT("logs", "eigenmt_prep_{condition}_chr{chr}_pc{pc}_{window}kb.log")
    shell:
        """
        module load r/{params.r_version}
        
        Rscript - <<'RSCRIPT' > {log} 2>&1
        
        library(data.table)
        library(dplyr)
        library(stringr)
        
        cat("Preparing EigenMT inputs\\n")
        cat("Condition: {wildcards.condition}\\n")
        cat("Chromosome: {wildcards.chr}\\n")
        cat("PC: {wildcards.pc}\\n")
        cat("Window: {wildcards.window}kb\\n")
        
        # Read adjusted RASQUAL results
        cat("Reading RASQUAL results from: {input.adjusted_results}\\n")
        rasqual <- fread("{input.adjusted_results}")
        
        cat("Total variants in adjusted results:", nrow(rasqual), "\\n")
        
        # Filter for current chromosome
        chrom_num <- gsub("chr", "", "{wildcards.chr}")
        rasqual_chr <- rasqual[Chromosome == "{wildcards.chr}" | Chromosome == chrom_num]
        
        cat("Variants for {wildcards.chr}:", nrow(rasqual_chr), "\\n")
        
        if (nrow(rasqual_chr) == 0) {{
            cat("WARNING: No variants found for chromosome {wildcards.chr}\\n")
            # Create empty output files
            empty_dt <- data.table()
            fwrite(empty_dt, "{output.gen_positions}", sep = "\\t")
            fwrite(empty_dt, "{output.genotypes}", sep = "\\t")
            fwrite(empty_dt, "{output.phe_positions}", sep = "\\t")
            fwrite(empty_dt, "{output.qtl}", sep = "\\t")
            quit(save = "no", status = 0)
        }}
        
        #############################
        # 1. gen.positions
        #############################
        cat("Creating gen.positions...\\n")
        gen_pos <- rasqual_chr %>%
            select(snp = rs_ID, chr_snp = Chromosome, pos = SNP_Position) %>%
            filter(snp != "SKIPPED", !is.na(snp), snp != "") %>%
            distinct(snp, .keep_all = TRUE)
        
        gen_pos$chr_snp <- str_remove(gen_pos$chr_snp, "^chr")
        
        fwrite(gen_pos, "{output.gen_positions}", sep = "\\t", 
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote gen.positions:", nrow(gen_pos), "variants\\n")
        
        #############################
        # 2. genotypes
        #############################
        cat("Reading genotype file: {input.traw}\\n")
        geno_dt <- fread("{input.traw}")
        
        cat("Genotype file dimensions:", nrow(geno_dt), "x", ncol(geno_dt), "\\n")
        
        # Extract SNP column and sample columns (columns 7 onward)
        genotypes <- geno_dt %>%
            select(SNP, 7:ncol(geno_dt)) %>%
            filter(SNP %in% gen_pos$snp)
        
        fwrite(genotypes, "{output.genotypes}", sep = "\\t",
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote genotypes:", nrow(genotypes), "variants x", ncol(genotypes)-1, "samples\\n")
        
        #############################
        # 3. phe.positions
        #############################
        cat("Creating phe.positions...\\n")
        phe_pos <- rasqual_chr %>%
            select(peak_id = Feature, chrom_probe = Chromosome) %>%
            distinct() %>%
            mutate(
                chrom_probe = str_remove(chrom_probe, "^chr"),
                tmp = str_remove(peak_id, paste0("^", "{wildcards.chr}", "_")),
                s1 = str_split_fixed(tmp, "_", 2)[, 1],
                s2 = str_split_fixed(tmp, "_", 2)[, 2]
            ) %>%
            select(peak_id, chrom_probe, s1, s2)
        
        fwrite(phe_pos, "{output.phe_positions}", sep = "\\t",
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote phe.positions:", nrow(phe_pos), "peaks\\n")
        
        #############################
        # 4. qtl
        #############################
        cat("Creating qtl file...\\n")
        qtl <- rasqual_chr %>%
            select(snp = rs_ID, peak = Feature, "p-value" = PValue) %>%
            filter(snp != "SKIPPED", !is.na(snp), snp != "", snp %in% gen_pos$snp)
        
        fwrite(qtl, "{output.qtl}", sep = "\\t",
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote qtl:", nrow(qtl), "associations\\n")
        
        cat("\\nEigenMT input preparation complete for {wildcards.chr}\\n")
        cat("Summary:\\n")
        cat("  - Variants:", nrow(gen_pos), "\\n")
        cat("  - Peaks:", nrow(phe_pos), "\\n")
        cat("  - Associations:", nrow(qtl), "\\n")
        
RSCRIPT
        """

rule run_eigenMT:
    input:
        qtl=rules.prepare_eigenmt_inputs.output.qtl,
        genotypes=rules.prepare_eigenmt_inputs.output.genotypes,
        gen_positions=rules.prepare_eigenmt_inputs.output.gen_positions,
        phe_positions=rules.prepare_eigenmt_inputs.output.phe_positions
    output:
        results=OUT("caQTL", "eigenMT", "results", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "eigenMT_chr{chr}_{condition}_pc{pc}.txt")
    params:
        eigenmt_dir=config['eigenmt_dir'],
        python_version=config['python'],
        window_snps=config['eigenmt_window_snps'],
        cis_dist=config['eigenmt_cis_dist'],
        var_thresh=config['eigenmt_var_thresh']
    threads: 1
    resources:
        mem_mb=10000,
        time="1-00:00:00"
    log:
        OUT("logs", "eigenmt_run_{condition}_chr{chr}_pc{pc}_{window}kb.log")
    shell:
        """
        module load python/{params.python_version}
        mkdir -p $(dirname {output.results})

        # Check if input files have data
        qtl_lines=$(wc -l < {input.qtl})

        if [ "$qtl_lines" -le 1 ]; then
            echo "WARNING: No data in QTL file for {wildcards.chr}" > {log}
            # Create empty output file
            echo -e "GENE\\tSNP\\tSTATISTIC\\tBETA\\tP_value" > {output.results}
        else
            python {params.eigenmt_dir}/eigenMT.py \
                --CHROM {wildcards.chr} \
                --QTL {input.qtl} \
                --GEN {input.genotypes} \
                --GENPOS {input.gen_positions} \
                --PHEPOS {input.phe_positions} \
                --var_thresh {params.var_thresh} \
                --cis_dist {params.cis_dist} \
                --OUT {output.results} \
                --window {params.window_snps} \
                2> {log}
        fi
        """
