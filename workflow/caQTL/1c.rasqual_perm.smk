#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts])

# Load analysis parameters from config
CONDITIONS = config['conditions']
VCF_GROUPS = config['vcf_groups']
CHROMOSOMES = config['chromosomes']
WINDOW_KB = config['window']
N_PERMUTATIONS = config['n_perm']

# Best PC for each condition
BEST_PC = config['best_pc']

# Get PC values to use (only best PC for each condition)
PC_VALUES = [BEST_PC[cond] for cond in CONDITIONS]

# Create mappings from config
CONDITION_TO_GROUP = dict(zip(CONDITIONS, VCF_GROUPS))

onsuccess:
    print("RASQUAL permutation analysis completed successfully!")
    print(f"Results available in: {OUT('caQTL', 'rasqual_permutation')}")

rule all:
    input:
        # Permutation RASQUAL results for all conditions, chromosomes (only best PC)
        expand(OUT("caQTL", "rasqual_permutation", "window_{window}kb", "pc{pc}",
                   "{condition}_{chr}_pc{pc}_{window}kb_perm{perm}.txt"),
               condition=CONDITIONS, chr=CHROMOSOMES,
               pc=PC_VALUES, window=WINDOW_KB,
               perm=range(1, N_PERMUTATIONS + 1)),
        #Filter
        expand(OUT("caQTL", "rasqual_permutation", "filtered_window_{window}kb", "pc{pc}",
            "{condition}_{chr}_pc{pc}_{window}kb_perm{perm}_filtered.txt"),condition=CONDITIONS, chr=CHROMOSOMES, 
            pc=PC_VALUES, window=WINDOW_KB, perm=range(1, N_PERMUTATIONS + 1)),
        # Lead SNPs
        expand(OUT("caQTL", "rasqual_permutation", "combined_window_{window}kb",
                   "{condition}_pc{pc}_{window}kb_perm{perm}_lead_snps.txt"),
               condition=CONDITIONS, pc=PC_VALUES, window=WINDOW_KB,
               perm=range(1, N_PERMUTATIONS + 1))


rule run_rasqual_permutation:
    input:
        peak_info=OUT("caQTL", "rasqual_input", "peakinfo", "peak_info_{chr}.txt"),
        expression=OUT("caQTL", "rasqual_input", "expression", "{chr}_{condition}.expression.bin"),
        size_factors=OUT("caQTL", "rasqual_input", "sizefactor", "{chr}_{condition}.size_factors.bin"),
        covariates=OUT("caQTL", "rasqual_input", "covar", "{condition}_covariates_pc{pc}.bin"),
        asvcf=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "03_asvcf",
                            f"{wc.chr}.{CONDITION_TO_GROUP[wc.condition]}.ASCounts.vcf.gz"),
        asvcf_tbi=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "03_asvcf",
                                f"{wc.chr}.{CONDITION_TO_GROUP[wc.condition]}.ASCounts.vcf.gz.tbi")
    output:
        result=OUT("caQTL", "rasqual_permutation", "window_{window}kb", "pc{pc}",
                   "{condition}_{chr}_pc{pc}_{window}kb_perm{perm}.txt")
    params:
        rasqual_bin=config['rasqual_dir'],
        tabix_ver=config['tabix_ver'],
        gsl_ver=config['gsl_ver'],
        vcf_dir=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "03_asvcf"),
        vcf_file=lambda wc: f"{wc.chr}.{CONDITION_TO_GROUP[wc.condition]}.ASCounts.vcf.gz",
        window_size=lambda wc: int(wc.window) * 1000,
        n_threads=5,
        maf_threshold=0.0,
        alpha=0.00,
        lambda_param=0.00,
        min_coverage=0.0
    threads: 10
    resources:
        mem_mb=5000,
        time="5-00:00:00"
    log:
        OUT("logs", "rasqual_perm_{condition}_{chr}_pc{pc}_{window}kb_perm{perm}.log")
    shell:
        """
        module load tabix/{params.tabix_ver}
        module load gsl/{params.gsl_ver}
        
        mkdir -p $(dirname {output.result})
        
        # Get sample number
        expression_txt=$(echo {input.expression} | sed 's/\\.bin$/.txt/')
        if [[ -f "$expression_txt" ]]; then
            sample_number=$(awk 'NR==1 {{print NF-1}}' "$expression_txt")
        else
            echo "Warning: Could not find $expression_txt, using default" >> {log}
            sample_number=21
        fi
        
        # Calculate total covariates (PC + 5 fixed: sex, batch, genoPCs 1-3)
        total_covariates=$(({wildcards.pc} + 5))
        
        echo "=== RASQUAL PERMUTATION ===" > {log}
        echo "Condition: {wildcards.condition}" >> {log}
        echo "Chromosome: {wildcards.chr}" >> {log}
        echo "PC: {wildcards.pc}" >> {log}
        echo "Permutation: {wildcards.perm}" >> {log}
        echo "Window: {params.window_size} bp" >> {log}
        echo "Samples: $sample_number" >> {log}
        echo "Covariates: $total_covariates" >> {log}
        echo "=========================" >> {log}
        
        # Initialize output with header
        echo -e "PeakID\\tChr\\tPeakStart\\tPeakEnd\\tPeakCenter\\tWindowStart\\tWindowEnd\\tNumCisSNPs\\tNumPeakSNPs\\tRASQUAL_Output" > {output.result}
        
        # Process each peak
        peak_rowNum=1
        while read peak; do
            # Skip header
            if [[ "$peak" == *"PeakID"* ]]; then
                continue
            fi
            
            # Parse peak info
            read -ra params <<< "$peak"
            peak_ID="${{params[0]}}"
            peak_chr="${{params[1]}}"
            peak_start="${{params[2]}}"
            peak_end="${{params[3]}}"
            
            # Calculate peak center
            peak_center=$(( (peak_start + peak_end) / 2 ))
            
            
            # Set window around center
            window_size={params.window_size}
            window_start=$(( peak_center - window_size ))
            (( window_start < 1 )) && window_start=1
            window_end=$(( peak_center + window_size ))
            
            # Calculate actual window size for RASQUAL -w parameter
            actual_window_size=$(( window_end - window_start ))
            
            # Count SNPs in regions
            num_cisSNPs=$(tabix {params.vcf_dir}/{params.vcf_file} ${{peak_chr}}:${{window_start}}-${{window_end}} | wc -l)
            num_peakSNPs=$(tabix {params.vcf_dir}/{params.vcf_file} ${{peak_chr}}:${{peak_start}}-${{peak_end}} | wc -l)
            
            echo "Peak #${{peak_rowNum}}: ${{peak_ID}} (perm {wildcards.perm})" >> {log}
            echo "  Window: ${{window_start}}-${{window_end}}, WindowSize: ${{actual_window_size}}" >> {log}
            echo "  cisSNPs: ${{num_cisSNPs}}, peakSNPs: ${{num_peakSNPs}}" >> {log}
            
            # Run RASQUAL with permutation (-r --force flags)
            rasqual_output=$(tabix {params.vcf_dir}/{params.vcf_file} ${{peak_chr}}:${{window_start}}-${{window_end}} | \\
                {params.rasqual_bin} \\
                    -s ${{peak_start}} -e ${{peak_end}} -j ${{peak_rowNum}} \\
                    -y {input.expression} -k {input.size_factors} \\
                    -n ${{sample_number}} -l ${{num_cisSNPs}} -m ${{num_peakSNPs}} \\
                    -f ${{peak_ID}} -x {input.covariates} -p ${{total_covariates}} \\
                    -w ${{actual_window_size}} -q {params.maf_threshold} \\
                    -a {params.alpha} -h {params.lambda_param} \\
                    --min-coverage-depth {params.min_coverage} \\
                    -r --force -z --n-threads {params.n_threads} 2>>{log})

            # Append to output
            echo -e "${{peak_ID}}\\t${{peak_chr}}\\t${{peak_start}}\\t${{peak_end}}\\t${{peak_center}}\\t${{window_start}}\\t${{window_end}}\\t${{num_cisSNPs}}\\t${{num_peakSNPs}}\\t${{rasqual_output}}" >> {output.result}
            
            peak_rowNum=$((peak_rowNum + 1))
            
        done < {input.peak_info}
        
        echo "Completed: {wildcards.chr} permutation {wildcards.perm}" >> {log}
        """

rule filter_permutation_results:
    input:
        raw_result=rules.run_rasqual_permutation.output.result
    output:
        filtered=OUT("caQTL", "rasqual_permutation", "filtered_window_{window}kb", "pc{pc}",
                    "{condition}_{chr}_pc{pc}_{window}kb_perm{perm}_filtered.txt")
    params:
        rsnp_threshold=config['rasqual_filters']['rsnp_correlation'],
        imputation_quality=config['rasqual_filters']['imputation_quality'],
        r_version=config['r']
    threads: 2
    resources:
        mem_mb=10000,
        time="1-00:00:00"
    log:
        OUT("logs", "filter_perm_{condition}_{chr}_pc{pc}_{window}kb_perm{perm}.log")
    shell:
        """
        awk -F'\\t' 'NF==25' {input.raw_result} > {output.filtered}.temp
        
        module load r/{params.r_version}
        Rscript - <<'RSCRIPT' > {log} 2>&1
        library(dplyr)
        library(data.table)

        cat("Starting RASQUAL filtering\\n")
        cat("Input: {output.filtered}.temp\\n")

        rasqual_col_names <- c(
          "Feature","rs_ID","Chromosome","SNP_Position",
          "Ref_Allele","Alt_allele","Allele_Frequency","HWE_Chi_Square",
          "Imputation_Quality","Log10_BH_Qvalue","Chi_square_statistic",
          "Effect_Size","Error_Rate","Ref_Bias","Overdispersion",
          "SNP_ID_Region","Feature_SNPs","Tested_SNPs","Null_Iterations",
          "Alt_Iterations","Random_Ties","Null_LogLikelihood",
          "Convergence","FSNP_Correlation","RSNP_Correlation"
        )

        # Read filtered file
        raw <- fread("{output.filtered}.temp", header=FALSE)
        cat("Total lines read:", nrow(raw), "\\n")

        if (nrow(raw) == 0) {{
          cat("Empty file, writing placeholder\\n")
          empty_dt <- data.frame(matrix(ncol=27, nrow=0))
          colnames(empty_dt) <- c(rasqual_col_names, "PValue", "Distance_From_Peak")
          fwrite(empty_dt, "{output.filtered}", sep="\\t")
          quit(save="no", status=0)
        }}

        rasqual_data <- raw
        colnames(rasqual_data) <- rasqual_col_names

        # Remove invalid positions
        rasqual_data <- rasqual_data %>% filter(SNP_Position != -1)

        # Cast numeric columns
        num_cols <- c("SNP_Position","Allele_Frequency","HWE_Chi_Square",
                      "Imputation_Quality","Log10_BH_Qvalue","Chi_square_statistic",
                      "Effect_Size","Error_Rate","Ref_Bias","Overdispersion",
                      "Feature_SNPs","Tested_SNPs","Null_Iterations","Alt_Iterations",
                      "Random_Ties","Null_LogLikelihood","Convergence",
                      "FSNP_Correlation","RSNP_Correlation")
        rasqual_data[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

        cat("Parsed rows:", nrow(rasqual_data), "\\n")

        # Filter by imputation quality
        rasqual_data <- rasqual_data %>% filter(Imputation_Quality >= {params.imputation_quality})
        cat("After imputation quality filter:", nrow(rasqual_data), "\\n")

        # Filter by RSNP correlation
        rasqual_data <- rasqual_data %>% filter(RSNP_Correlation >= {params.rsnp_threshold})
        cat("After RSNP correlation filter:", nrow(rasqual_data), "\\n")

        # Calculate PValue
        rasqual_data$PValue <- pchisq(rasqual_data$Chi_square_statistic, df=1, lower.tail=FALSE)
        rasqual_data$PValue[is.na(rasqual_data$PValue) | rasqual_data$PValue < 0] <- 1

        # Calculate Distance_From_Peak using data.table syntax
        rasqual_data[, c("chr_temp", "Peak_Start", "Peak_End") := tstrsplit(Feature, "_", fixed=TRUE)]
        rasqual_data[, Peak_Start := as.numeric(Peak_Start)]
        rasqual_data[, Peak_End := as.numeric(Peak_End)]
        rasqual_data[, Peak_Center := round((Peak_End - Peak_Start)/2 + Peak_Start)]
        rasqual_data[, Distance_From_Peak := abs(SNP_Position - Peak_Center)]
        rasqual_data[, c("chr_temp", "Peak_Start", "Peak_End", "Peak_Center") := NULL]

        # Save out
        fwrite(rasqual_data, "{output.filtered}", sep="\\t", quote=FALSE)

        cat("Filtered results written. Final rows:", nrow(rasqual_data), "\\n")
RSCRIPT
        
        # Clean up temp file
        rm -f {output.filtered}.temp
        """

rule combine_permutation_results:
    input:
        filtered_files=expand(OUT("caQTL", "rasqual_permutation", "filtered_window_{{window}}kb", "pc{{pc}}",
                              "{{condition}}_{chr}_pc{{pc}}_{{window}}kb_perm{{perm}}_filtered.txt"),
                          chr=CHROMOSOMES)
    output:
        combined=OUT("caQTL", "rasqual_permutation", "combined_window_{window}kb",
                    "{condition}_pc{pc}_{window}kb_perm{perm}_combined.txt"),
        stats=OUT("caQTL", "rasqual_permutation", "combined_window_{window}kb",
                "{condition}_pc{pc}_{window}kb_perm{perm}_stats.txt")
    params:
        r_version=config['r']
    threads: 4
    resources:
        mem_mb=50000,
        time="1-00:00:00"
    log:
        OUT("logs", "combine_perm_{condition}_pc{pc}_{window}kb_perm{perm}.log")
    shell:
        """
        module load r/{params.r_version}

        # Create file list from Snakemake input
        echo "{input.filtered_files}" | tr ' ' '\\n' > {log}.files.txt

        cat > {log}.R <<'RSCRIPT'
library(data.table)
library(dplyr)

cat("Starting chromosome combination\\n")

input_files <- readLines("{log}.files.txt")
cat("Number of input files:", length(input_files), "\\n")

combined_data <- NULL
stats_list <- list()

for (i in seq_along(input_files)) {{
    file_path <- input_files[i]
    chr_name <- basename(file_path)

    cat("Processing:", chr_name, "\\n")

    tryCatch({{
        chr_data <- fread(file_path, header = TRUE, stringsAsFactors = FALSE)

        if (nrow(chr_data) > 0) {{
            stats_list[[chr_name]] <- data.frame(
                chromosome = chr_name,
                n_variants = nrow(chr_data),
                n_peaks = length(unique(chr_data$Feature))
            )

            if (is.null(combined_data)) {{
                combined_data <- chr_data
            }} else {{
                combined_data <- rbindlist(list(combined_data, chr_data),
                                          use.names = TRUE, fill = TRUE)
            }}

            cat("  Rows:", nrow(chr_data), "\\n")
        }}
    }}, error = function(e) {{
        cat("  Error:", e$message, "\\n")
    }})
}}

if (is.null(combined_data) || nrow(combined_data) == 0) {{
    cat("No data to combine\\n")
    empty_dt <- data.table(
        Feature = character(), rs_ID = character(), Chromosome = character(),
        SNP_Position = numeric(), Ref_Allele = character(), Alt_allele = character(),
        Allele_Frequency = numeric(), HWE_Chi_Square = numeric(),
        Imputation_Quality = numeric(), Log10_BH_Qvalue = numeric(),
        Chi_square_statistic = numeric(), Effect_Size = numeric(),
        Error_Rate = numeric(), Ref_Bias = numeric(), Overdispersion = numeric(),
        SNP_ID_Region = character(), Feature_SNPs = integer(), Tested_SNPs = integer(),
        Null_Iterations = integer(), Alt_Iterations = integer(),
        Random_Ties = numeric(), Null_LogLikelihood = numeric(),
        Convergence = integer(), FSNP_Correlation = numeric(),
        RSNP_Correlation = numeric(), PValue = numeric(),
        Distance_From_Peak = numeric()
    )
    fwrite(empty_dt, "{output.combined}", sep = "\\t")

    empty_stats <- data.table(
        chromosome = character(),
        n_variants = integer(),
        n_peaks = integer()
    )
    fwrite(empty_stats, "{output.stats}", sep = "\\t")
}} else {{
    fwrite(combined_data, "{output.combined}", sep = "\\t",
           col.names = TRUE, row.names = FALSE, quote = FALSE)

    if (length(stats_list) > 0) {{
        stats_df <- rbindlist(stats_list)
        fwrite(stats_df, "{output.stats}", sep = "\\t",
               col.names = TRUE, row.names = FALSE, quote = FALSE)
    }}

    cat("\\nCombination complete\\n")
    cat("Total variants:", nrow(combined_data), "\\n")
    cat("Total unique peaks:", length(unique(combined_data$Feature)), "\\n")
}}
RSCRIPT

        Rscript {log}.R 2>&1 | tee {log}
        rm -f {log}.files.txt {log}.R
        """


rule adjust_tied_perm:
    input:
        combined=rules.combine_permutation_results.output.combined
    output:
        adjusted=OUT("caQTL", "rasqual_permutation", "combined_window_{window}kb",
                "{condition}_pc{pc}_{window}kb_perm{perm}_adjusted.txt")
    params:
        r_version=config['r']
    threads: 4
    resources:
        mem_mb=75000,
        time="5-00:00:00"
    log:
        OUT("logs", "adjusted_perm_{condition}_pc{pc}_{window}kb_perm{perm}.log")
    shell:
        """
        module load r/{params.r_version}

        Rscript - <<'RSCRIPT' > {log} 2>&1

        library(data.table)
        library(dplyr)

        cat("Starting tied p-value adjustment\\n")

        combined_data <- fread("{input.combined}", header = TRUE, stringsAsFactors = FALSE)

        if (nrow(combined_data) == 0) {{
            cat("No data to adjust\\n")
            file.copy("{input.combined}", "{output.adjusted}", overwrite = TRUE)
            quit(save = "no", status = 0)
        }}

        cat("Input rows:", nrow(combined_data), "\\n")
        cat("Unique peaks:", length(unique(combined_data$Feature)), "\\n")

        features <- unique(combined_data$Feature)
        cat("Processing", length(features), "features for tie adjustment\\n")

        adjusted_data <- copy(combined_data)
        features_with_ties <- 0
        total_ties_adjusted <- 0

        for (j in seq_along(features)) {{
            if (j %% 1000 == 0 || j == length(features)) {{
                cat("  Processed", j, "of", length(features), "features\\n")
            }}

            current_feature <- features[j]
            feature_data <- adjusted_data[Feature == current_feature]

            if (nrow(feature_data) <= 1) next

            setorder(feature_data, PValue, Distance_From_Peak)

            dup_indices <- which(duplicated(feature_data$PValue))

            if (length(dup_indices) > 0) {{
                feature_data$PValue[dup_indices] <- feature_data$PValue[dup_indices] + 1e-15
                features_with_ties <- features_with_ties + 1
                total_ties_adjusted <- total_ties_adjusted + length(dup_indices)

                adjusted_data[Feature == current_feature, PValue := feature_data$PValue]
            }}
        }}

        cat("Features with ties:", features_with_ties, "out of", length(features),
            sprintf("(%.2f%%)\\n", 100 * features_with_ties / length(features)))
        cat("Total ties adjusted:", total_ties_adjusted, "\\n")

        fwrite(adjusted_data, "{output.adjusted}", sep = "\\t",
               col.names = TRUE, row.names = FALSE, quote = FALSE)

        cat("Adjusted results written\\n")

RSCRIPT
        """


rule cleanup_matrix:
    input:
        adjusted=rules.adjust_tied_perm.output.adjusted
    output:
        lead_snps=OUT("caQTL", "rasqual_permutation", "combined_window_{window}kb",
                "{condition}_pc{pc}_{window}kb_perm{perm}_lead_snps.txt")
    params:
        r_version=config['r']
    threads: 2
    resources:
        mem_mb=50000,
        time="2-00:00:00"
    log:
        OUT("logs", "lead_snps_perm_{condition}_pc{pc}_{window}kb_perm{perm}.log")
    shell:
        """
        module load r/{params.r_version}

        Rscript - <<'RSCRIPT' > {log} 2>&1

        library(data.table)

        cat("Preparing ALL adjusted SNPs for eigenMT (Wnt approach)\\n")

        adjusted_data <- fread("{input.adjusted}", header = TRUE, stringsAsFactors = FALSE)

        if (nrow(adjusted_data) == 0) {{
            cat("No data to process\\n")
            empty_dt <- data.table(
                Feature = character(), rs_ID = character(), Chromosome = character(),
                SNP_Position = numeric(), Ref_Allele = character(), Alt_allele = character(),
                Allele_Frequency = numeric(), HWE_Chi_Square = numeric(),
                Imputation_Quality = numeric(), Log10_BH_Qvalue = numeric(),
                Chi_square_statistic = numeric(), Effect_Size = numeric(),
                Error_Rate = numeric(), Ref_Bias = numeric(), Overdispersion = numeric(),
                SNP_ID_Region = character(), Feature_SNPs = integer(), Tested_SNPs = integer(),
                Null_Iterations = integer(), Alt_Iterations = integer(),
                Random_Ties = numeric(), Null_LogLikelihood = numeric(),
                Convergence = integer(), FSNP_Correlation = numeric(),
                RSNP_Correlation = numeric(), PValue = numeric(),
                Distance_From_Peak = numeric()
            )
            fwrite(empty_dt, "{output.lead_snps}", sep = "\\t")
            quit(save = "no", status = 0)
        }}

        cat("Input rows:", nrow(adjusted_data), "\\n")
        cat("Unique peaks:", length(unique(adjusted_data$Feature)), "\\n")

        fwrite(adjusted_data, "{output.lead_snps}", sep = "\\t",
               col.names = TRUE, row.names = FALSE, quote = FALSE)

        cat("All adjusted SNPs written for eigenMT\\n")

RSCRIPT
        """

