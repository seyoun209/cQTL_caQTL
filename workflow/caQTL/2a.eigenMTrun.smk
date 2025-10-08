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
                   "gen.positions_{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        expand(OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "genotypes_{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        expand(OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "phe.positions_{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        expand(OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                   "{condition}", "pc{pc}",
                   "qtl_{chr}_{condition}_pc{pc}.txt"),
               condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB),
        #eigenMT results
        expand(OUT("caQTL", "eigenMT", "results", "window_{window}kb","{condition}", "pc{pc}",
                   "eigenMT_{chr}_{condition}_pc{pc}.txt"),condition=CONDITIONS, pc=PC_VALUES, chr=CHROMOSOMES, window=WINDOW_KB)

rule prepare_eigenmt_inputs:
    input:
        adjusted_results=OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                           "{condition}_pc{pc}_{window}kb_adjusted.txt"),
        traw=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "04_recode",
                           "{chr}." + CONDITION_TO_GROUP[wc.condition] + ".traw")
    output:
        gen_positions=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                         "{condition}", "pc{pc}",
                         "gen.positions_{chr}_{condition}_pc{pc}.txt"),
        genotypes=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                     "{condition}", "pc{pc}",
                     "genotypes_{chr}_{condition}_pc{pc}.txt"),
        phe_positions=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
                         "{condition}", "pc{pc}",
                         "phe.positions_{chr}_{condition}_pc{pc}.txt"),
        qtl=OUT("caQTL", "eigenMT", "inputs", "window_{window}kb",
               "{condition}", "pc{pc}",
               "qtl_{chr}_{condition}_pc{pc}.txt")
    params:
        r_version=config['r']
    threads: 2
    resources:
        mem_mb=55000,
        time="1:00:00"
    log:
        OUT("logs", "eigenmt_prep_{condition}_{chr}_pc{pc}_{window}kb.log")
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
        
        #############################
        # Load RASQUAL results
        #############################
        cat("Reading RASQUAL results from: {input.adjusted_results}\\n")
        RASQUAL <- fread("{input.adjusted_results}")
        
        # Filter for current chromosome
        chrom_num <- gsub("chr", "", "{wildcards.chr}")
        RASQUAL_chr <- RASQUAL[Chromosome == "{wildcards.chr}" | Chromosome == chrom_num]
        
        cat("Total variants for {wildcards.chr}:", nrow(RASQUAL_chr), "\\n")
        
        if (nrow(RASQUAL_chr) == 0) {{
            cat("WARNING: No variants found for chromosome {wildcards.chr}\\n")
            stop("No data for this chromosome")
        }}
        
        #############################
        # 1. gen.positions
        #############################
        cat("Creating gen.positions...\\n")
        Gen.pos <- RASQUAL_chr %>%
            select(snp = rs_ID, chr_snp = Chromosome, pos = SNP_Position) %>%
            filter(snp != "SKIPPED") %>%
            distinct(snp, .keep_all = TRUE)
        
        Gen.pos$chr_snp <- str_remove(Gen.pos$chr_snp, "^chr")
        
        fwrite(Gen.pos, "{output.gen_positions}", sep = "\\t", 
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote gen.positions:", nrow(Gen.pos), "variants\\n")
        
        #############################
        # 2. genotypes
        #############################
        cat("Reading genotype file: {input.traw}\\n")
        geno_dt <- fread("{input.traw}")
        
        num <- ncol(geno_dt)
        genotypes <- geno_dt[, c(2, 7:num), with = FALSE]
        setnames(genotypes, old = colnames(genotypes)[1], new = "SNP")
        genotypes <- genotypes[SNP %in% Gen.pos$snp]
        
        fwrite(genotypes, "{output.genotypes}", sep = "\\t",
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote genotypes:", nrow(genotypes), "variants x", ncol(genotypes)-1, "samples\\n")
        
        #############################
        # 3. phe.positions
        #############################
        cat("Creating phe.positions...\\n")
        Phe.pos <- RASQUAL_chr %>%
            select(peak_id = Feature, chrom_probe = Chromosome) %>%
            distinct()
        
        # Remove chromosome prefix from peak names
        chrop <- paste0("{wildcards.chr}", "_")
        A <- gsub(chrop, "", Phe.pos$peak_id)
        A <- str_split_fixed(A, "_", 2)
        
        Phe.pos <- cbind(Phe.pos, A)
        colnames(Phe.pos) <- c("peak_id", "chrom_probe", "s1", "s2")
        Phe.pos$chrom_probe <- str_remove(Phe.pos$chrom_probe, "^chr")
        
        fwrite(Phe.pos, "{output.phe_positions}", sep = "\\t",
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote phe.positions:", nrow(Phe.pos), "peaks\\n")
        
        #############################
        # 4. qtl
        #############################
        cat("Creating qtl file...\\n")
        qtl <- RASQUAL_chr %>%
            select(snp = rs_ID, peak = Feature, "p-value" = PValue) %>%
            filter(snp != "SKIPPED")
        
        fwrite(qtl, "{output.qtl}", sep = "\\t",
               quote = FALSE, row.names = FALSE, col.names = TRUE)
        cat("Wrote qtl:", nrow(qtl), "associations\\n")
        
        cat("\\nEigenMT input preparation complete for {wildcards.chr}\\n")
        cat("Summary:\\n")
        cat("  - Variants:", nrow(Gen.pos), "\\n")
        cat("  - Peaks:", nrow(Phe.pos), "\\n")
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
                   "eigenMT_{chr}_{condition}_pc{pc}.txt")
    params:
        eigenmt_dir=config['eigenmt_dir'],
        python_version=config['python'],
        window_snps=config['eigenmt_window_snps'],
        cis_dist=config['eigenmt_cis_dist'],
        var_thresh=config['eigenmt_var_thresh']
    threads: 1
    resources:
        mem_mb=70000,
        time="1-00:00:00"
    log:
        OUT("logs", "eigenmt_run_{condition}_{chr}_pc{pc}_{window}kb.log")
    shell:
        """
        module load python/{params.python_version}
        mkdir -p $(dirname {output.results})

        # Remove 'chr' prefix from chromosome for eigenMT
        CHROM=$(echo {wildcards.chr} | sed 's/chr//')

        # Check if input files have data
        qtl_lines=$(wc -l < {input.qtl})

        if [ "$qtl_lines" -le 1 ]; then
            echo "WARNING: No data in QTL file for {wildcards.chr}" > {log}
            # Create empty output file
            echo -e "GENE\\tSNP\\tSTATISTIC\\tBETA\\tP_value" > {output.results}
        else
            python {params.eigenmt_dir}/eigenMT.py \
                --CHROM "$CHROM" \
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
