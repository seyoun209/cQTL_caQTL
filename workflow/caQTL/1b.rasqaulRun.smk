#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RASQUAL Execution Pipeline - 1b.rasqualRun.smk
Clean version for caQTL mapping
"""

import os
import pandas as pd

def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

# Load analysis parameters from config
CONDITIONS = config['conditions']
VCF_GROUPS = config['vcf_groups']
CHROMOSOMES = config['chromosomes']
PC_VALUES = config['atac_pcs']
WINDOW_KB = config['window']

# Create mappings from config
CONDITION_TO_GROUP = dict(zip(CONDITIONS, VCF_GROUPS))
GROUP_TO_CONDITION = dict(zip(VCF_GROUPS, CONDITIONS))

onsuccess:
    print("RASQUAL analysis completed successfully!")
    print(f"Results available in: {OUT('caQTL', 'rasqual_output')}")

rule all:
    input:
        # RASQUAL results for all conditions, chromosomes, and PC values
        expand(OUT("caQTL", "rasqual_output", "window_{window}kb", "pc{pc}",
                   "{condition}_{chr}_pc{pc}_{window}kb.txt"),condition=CONDITIONS, chr=CHROMOSOMES, pc=PC_VALUES, window=WINDOW_KB),

        # Combined results (all chromosomes merged)
        expand(OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                   "{condition}_pc{pc}_{window}kb_combined.txt"),condition=CONDITIONS, pc=PC_VALUES, window=WINDOW_KB),
        
        # Adjusted results (tied p-values fixed)
        expand(OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                   "{condition}_pc{pc}_{window}kb_adjusted.txt"),condition=CONDITIONS, pc=PC_VALUES, window=WINDOW_KB),
        
        # eigenMT input snps
        expand(OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                   "{condition}_pc{pc}_{window}kb_lead_snps.txt"),condition=CONDITIONS, pc=PC_VALUES, window=WINDOW_KB)

rule run_rasqual:
    input:
        # Peak info file
        peak_info=OUT("caQTL", "rasqual_input", "peakinfo", "peak_info_{chr}.txt"),
        # Expression counts (binary)
        expression=OUT("caQTL", "rasqual_input", "expression", "{chr}_{condition}.expression.bin"),
        # Size factors (binary)
        size_factors=OUT("caQTL", "rasqual_input", "sizefactor", "{chr}_{condition}.size_factors.bin"),
        # Covariates (binary)
        covariates=OUT("caQTL", "rasqual_input", "covar", "{condition}_covariates_pc{pc}.bin"),
        # AS-VCF file
        asvcf=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "03_asvcf",
                            f"{wc.chr}.{CONDITION_TO_GROUP[wc.condition]}.ASCounts.vcf.gz"),
        asvcf_tbi=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "03_asvcf",
                                f"{wc.chr}.{CONDITION_TO_GROUP[wc.condition]}.ASCounts.vcf.gz.tbi")
    output:
        result=OUT("caQTL", "rasqual_output", "window_{window}kb", "pc{pc}",
                   "{condition}_{chr}_pc{pc}_{window}kb.txt")
    params:
        rasqual_bin=config['rasqual_dir'],
        tabix_ver=config['tabix_ver'],
        gsl_ver=config['gsl_ver'],
        vcf_dir=lambda wc: OUT("geno", CONDITION_TO_GROUP[wc.condition], "03_asvcf"),
        vcf_file=lambda wc: f"{wc.chr}.{CONDITION_TO_GROUP[wc.condition]}.ASCounts.vcf.gz",
        window_size=lambda wc: int(wc.window) * 1000,  # Convert kb to bp
        n_threads=10,
        # RASQUAL parameters
        maf_threshold=0.0,
        alpha=0.00,
        lambda_param=0.00,
        min_coverage=0.0
    threads: 10
    resources:
        mem_mb=5000,
        time="3-00:00:00"
    log:
        OUT("logs", "rasqual_{condition}_{chr}_pc{pc}_{window}kb.log")
    shell:
        """
        module load tabix/{params.tabix_ver}
        module load gsl/{params.gsl_ver}
        
        # Create output directory
        mkdir -p $(dirname {output.result})
        
        # Get sample number from expression text file
        expression_txt=$(echo {input.expression} | sed 's/\\.bin$/.txt/')
        if [[ -f "$expression_txt" ]]; then
            sample_number=$(awk 'NR==1 {{print NF-1}}' "$expression_txt")
        else
            echo "Warning: Could not find $expression_txt, using default sample number"
            sample_number=21
        fi
        
        # Calculate total covariates (PC value + 5) sex, atac_submission_batch, genopc1-3
        total_covariates=$(({wildcards.pc} + 5))
        
        echo "=== RASQUAL Parameters ===" > {log}
        echo "Condition: {wildcards.condition}" >> {log}
        echo "Chromosome: {wildcards.chr}" >> {log}
        echo "PC value: {wildcards.pc}" >> {log}
        echo "Window size: {params.window_size} bp" >> {log}
        echo "Sample number: $sample_number" >> {log}
        echo "Total covariates: $total_covariates" >> {log}
        echo "=========================" >> {log}
        
        # Initialize output file with header
        echo -e "PeakID\\tChr\\tPeakStart\\tPeakEnd\\tPeakCenter\\tWindowStart\\tWindowEnd\\tNumCisSNPs\\tNumPeakSNPs\\tRASQUAL_Output" > {output.result}
        
        # Process each peak
        peak_rowNum=1
        while read peak; do
            # Skip header line
            if [[ "$peak" == *"PeakID"* ]]; then
                continue
            fi
            
            # Parse peak information
            read -ra params <<< "$peak"
            peak_ID="${{params[0]}}"
            peak_chr="${{params[1]}}"
            peak_start="${{params[2]}}"
            peak_end="${{params[3]}}"
            
            # Calculate peak center
            peak_center=$(( (peak_start + peak_end) / 2 ))
            
            # Set window around peak center
            
            window_size={params.window_size}
            window_start=$(( peak_center - window_size / 2 ))
            (( window_start < 1 )) && window_start=1
            window_end=$(( peak_center + window_size / 2 ))


            # Count SNPs in window and peak
            num_cisSNPs=$(tabix {params.vcf_dir}/{params.vcf_file} ${{peak_chr}}:${{window_start}}-${{window_end}} | wc -l)
            num_peakSNPs=$(tabix {params.vcf_dir}/{params.vcf_file} ${{peak_chr}}:${{peak_start}}-${{peak_end}} | wc -l)
            
            echo "Processing peak #${{peak_rowNum}}: ${{peak_ID}} (${{peak_chr}}:${{peak_start}}-${{peak_end}})" >> {log}
            echo "  Window: ${{window_start}}-${{window_end}}, SNPs: ${{num_cisSNPs}}" >> {log}
            
            # Run RASQUAL
            rasqual_output=$(tabix {params.vcf_dir}/{params.vcf_file} ${{peak_chr}}:${{window_start}}-${{window_end}} | \\
                {params.rasqual_bin} \\
                    -s ${{peak_start}} -e ${{peak_end}} -j ${{peak_rowNum}} \\
                    -y {input.expression} -k {input.size_factors} \\
                    -n ${{sample_number}} -l ${{num_cisSNPs}} -m ${{num_peakSNPs}} \\
                    -f ${{peak_ID}} -x {input.covariates} -p ${{total_covariates}} \\
                    -w {params.window_size} -q {params.maf_threshold} \\
                    -a {params.alpha} -h {params.lambda_param} \\
                    --min-coverage-depth {params.min_coverage} \\
                    -z --n-threads {params.n_threads} 2>>{log})
            
            # Append result to output file
            echo -e "${{peak_ID}}\\t${{peak_chr}}\\t${{peak_start}}\\t${{peak_end}}\\t${{peak_center}}\\t${{window_start}}\\t${{window_end}}\\t${{num_cisSNPs}}\\t${{num_peakSNPs}}\\t${{rasqual_output}}" >> {output.result}
            
            peak_rowNum=$((peak_rowNum + 1))
            
        done < {input.peak_info}
        
        echo "RASQUAL processing completed. Total peaks processed: $((peak_rowNum-1))" >> {log}
        """


rule filter_chromosome_results:
    input:
        raw_result=rules.run_rasqual.output.result
    output:
        filtered=OUT("caQTL", "rasqual_output", "filtered_window_{window}kb", "pc{pc}",
                    "{condition}_{chr}_pc{pc}_{window}kb_filtered.txt")
    params:
        rsnp_threshold=config['rasqual_filters']['rsnp_correlation'],
        imputation_quality=config['rasqual_filters']['imputation_quality'],
        r_version=config['r']
    threads: 2
    resources:
        mem_mb=8000,
        time="5:00:00"
    log:
        OUT("logs", "filter_{condition}_{chr}_pc{pc}_{window}kb.log")
    shell:
        """
        module load r/{params.r_version}
        
        Rscript - <<'RSCRIPT' > {log} 2>&1
        
        library(data.table)
        library(dplyr)
        library(stringr)
        
        cat("Starting chromosome filtering\\n")
        cat("Input: {input.raw_result}\\n")
        cat("RSNP threshold: {params.rsnp_threshold}\\n")
        cat("Imputation quality: {params.imputation_quality}\\n")
        
        # Read raw RASQUAL output - read as text first to handle variable columns
        rasqual_lines <- readLines("{input.raw_result}")
        
        if (length(rasqual_lines) <= 1) {{
            cat("No data to process (only header or empty)\\n")
            # Create empty output with proper headers
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
            fwrite(empty_dt, "{output.filtered}", sep = "\\t")
            quit(save = "no", status = 0)
        }}
        
        # Parse each line manually
        header <- strsplit(rasqual_lines[1], "\\t")[[1]]
        cat("Header columns:", length(header), "\\n")
        cat("Header:", paste(header, collapse = ", "), "\\n")
        
        # Process data lines
        data_lines <- rasqual_lines[-1]  # Skip header
        cat("Initial data rows:", length(data_lines), "\\n")
        
        if (length(data_lines) == 0) {{
            cat("No data rows to process\\n")
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
            fwrite(empty_dt, "{output.filtered}", sep = "\\t")
            quit(save = "no", status = 0)
        }}
        
        # Parse each data line: first 9 columns are metadata, rest is RASQUAL output
        parsed_data <- lapply(data_lines, function(line) {{
            parts <- strsplit(line, "\\t")[[1]]
            
            # First 9 columns are the metadata from our shell script
            metadata <- parts[1:9]
            
            # Everything after column 9 is the RASQUAL output (tab-delimited)
            rasqual_output <- parts[10:length(parts)]
            
            return(list(
                metadata = metadata,
                rasqual = rasqual_output
            ))
        }})
        
        # Extract metadata
        metadata_matrix <- do.call(rbind, lapply(parsed_data, function(x) x$metadata))
        colnames(metadata_matrix) <- header[1:9]
        metadata_df <- as.data.frame(metadata_matrix, stringsAsFactors = FALSE)
        
        # Extract RASQUAL columns
        rasqual_matrix <- do.call(rbind, lapply(parsed_data, function(x) {{
            r <- x$rasqual
            # Ensure all rows have same length
            length(r) <- 25  # RASQUAL outputs 25 columns
            return(r)
        }}))
        
        # Define RASQUAL column names
        rasqual_col_names <- c("Feature", "rs_ID", "Chromosome", "SNP_Position", "Ref_Allele",
                               "Alt_allele", "Allele_Frequency", "HWE_Chi_Square",
                               "Imputation_Quality", "Log10_BH_Qvalue",
                               "Chi_square_statistic", "Effect_Size",
                               "Error_Rate", "Ref_Bias",
                               "Overdispersion", "SNP_ID_Region", "Feature_SNPs",
                               "Tested_SNPs", "Null_Iterations", "Alt_Iterations",
                               "Random_Ties", "Null_LogLikelihood", "Convergence",
                               "FSNP_Correlation", "RSNP_Correlation")
        
        colnames(rasqual_matrix) <- rasqual_col_names
        rasqual_df <- as.data.frame(rasqual_matrix, stringsAsFactors = FALSE)
        
        # Convert to data.table
        rasqual_parsed <- as.data.table(rasqual_df)
        
        # Convert numeric columns
        numeric_cols <- c("SNP_Position", "Allele_Frequency", "HWE_Chi_Square",
                         "Imputation_Quality", "Log10_BH_Qvalue", "Chi_square_statistic",
                         "Effect_Size", "Error_Rate", "Ref_Bias", "Overdispersion",
                         "Feature_SNPs", "Tested_SNPs", "Null_Iterations",
                         "Alt_Iterations", "Random_Ties", "Null_LogLikelihood",
                         "Convergence", "FSNP_Correlation", "RSNP_Correlation")
        
        for (col in numeric_cols) {{
            if (col %in% names(rasqual_parsed)) {{
                rasqual_parsed[[col]] <- as.numeric(rasqual_parsed[[col]])
            }}
        }}
        
        cat("Parsed rows:", nrow(rasqual_parsed), "\\n")
        
        # Filter out invalid SNP positions
        rasqual_parsed <- rasqual_parsed[!is.na(SNP_Position) & SNP_Position != -1]
        cat("After removing invalid positions:", nrow(rasqual_parsed), "\\n")
        
        if (nrow(rasqual_parsed) == 0) {{
            cat("No valid data after filtering\\n")
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
            fwrite(empty_dt, "{output.filtered}", sep = "\\t")
            quit(save = "no", status = 0)
        }}
        
        # Apply RSNP correlation filter
        rasqual_parsed <- rasqual_parsed[RSNP_Correlation >= {params.rsnp_threshold}]
        cat("After RSNP filter (>= {params.rsnp_threshold}):", nrow(rasqual_parsed), "\\n")
        
        if (nrow(rasqual_parsed) == 0) {{
            cat("No data after RSNP filter\\n")
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
            fwrite(empty_dt, "{output.filtered}", sep = "\\t")
            quit(save = "no", status = 0)
        }}
        
        # Apply imputation quality filter
        rasqual_parsed <- rasqual_parsed[Imputation_Quality >= {params.imputation_quality}]
        cat("After imputation quality filter (>= {params.imputation_quality}):", nrow(rasqual_parsed), "\\n")
        
        if (nrow(rasqual_parsed) == 0) {{
            cat("No data after imputation quality filter\\n")
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
            fwrite(empty_dt, "{output.filtered}", sep = "\\t")
            quit(save = "no", status = 0)
        }}
        
        # Calculate p-value
        rasqual_parsed[, PValue := pchisq(Chi_square_statistic, df = 1, lower.tail = FALSE)]
        rasqual_parsed[, PValue := ifelse(is.na(PValue) | PValue < 0, 1, PValue)]
        
        # Calculate distance from peak center using the Feature column
        rasqual_parsed[, c("chr_temp", "Peak_Start", "Peak_End") := tstrsplit(Feature, "_", keep = 1:3)]
        rasqual_parsed[, `:=`(
            Peak_Start = as.numeric(Peak_Start),
            Peak_End = as.numeric(Peak_End)
        )]
        rasqual_parsed[, Peak_Center := round((Peak_End - Peak_Start)/2 + Peak_Start)]
        rasqual_parsed[, Distance_From_Peak := abs(SNP_Position - Peak_Center)]
        rasqual_parsed[, c("chr_temp", "Peak_Start", "Peak_End", "Peak_Center") := NULL]
        
        # Write filtered results
        fwrite(rasqual_parsed, "{output.filtered}", sep = "\\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
        
        cat("Filtered results written\\n")
        cat("Final rows:", nrow(rasqual_parsed), "\\n")
        
RSCRIPT
        """

rule combine_chromosomes:
    input:
        filtered_files=lambda wc: expand(
            OUT("caQTL", "rasqual_output", "filtered_window_{{window}}kb", "pc{{pc}}",
                "{{condition}}_{chr}_pc{{pc}}_{{window}}kb_filtered.txt"),
            chr=CHROMOSOMES)
    output:
        combined=OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                    "{condition}_pc{pc}_{window}kb_combined.txt"),
        stats=OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                 "{condition}_pc{pc}_{window}kb_stats.txt")
    params:
        r_version=config['r']
    threads: 4
    resources:
        mem_mb=32000,
        time="5:00:00"
    log:
        OUT("logs", "combine_{condition}_pc{pc}_{window}kb.log")
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

rule adjust_tied_pvalues:
    input:
        combined=rules.combine_chromosomes.output.combined
    output:
        adjusted=OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                    "{condition}_pc{pc}_{window}kb_adjusted.txt")
    params:
        r_version=config['r']
    threads: 4
    resources:
        mem_mb=32000,
        time="3:00:00"
    log:
        OUT("logs", "adjust_{condition}_pc{pc}_{window}kb.log")
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

rule identify_lead_snps:
    input:
        adjusted=rules.adjust_tied_pvalues.output.adjusted
    output:
        lead_snps=OUT("caQTL", "rasqual_output", "combined_window_{window}kb",
                     "{condition}_pc{pc}_{window}kb_lead_snps.txt")
    params:
        r_version=config['r']
    threads: 2
    resources:
        mem_mb=16000,
        time="3:00:00"
    log:
        OUT("logs", "lead_snps_{condition}_pc{pc}_{window}kb.log")
    shell:
        """
        module load r/{params.r_version}

        Rscript - <<'RSCRIPT' > {log} 2>&1

        library(data.table)
        library(dplyr)

        cat("Starting lead SNP identification\\n")

        adjusted_data <- fread("{input.adjusted}", header = TRUE, stringsAsFactors = FALSE)

        if (nrow(adjusted_data) == 0) {{
            cat("No data to process\\n")
            # Create empty output with headers
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

        features <- unique(adjusted_data$Feature)
        cat("Identifying lead SNPs for", length(features), "features\\n")

        lead_snps <- data.table()

        for (i in seq_along(features)) {{
            if (i %% 1000 == 0 || i == length(features)) {{
                cat("  Processed", i, "of", length(features), "features\\n")
            }}

            current_feature <- features[i]
            feature_data <- adjusted_data[Feature == current_feature]

            setorder(feature_data, PValue, Distance_From_Peak)

            lead_snp <- feature_data[1, ]
            lead_snps <- rbindlist(list(lead_snps, lead_snp), use.names = TRUE, fill = TRUE)
        }}

        fwrite(lead_snps, "{output.lead_snps}", sep = "\\t",
               col.names = TRUE, row.names = FALSE, quote = FALSE)

        cat("Lead SNPs identified:", nrow(lead_snps), "\\n")
        cat("Lead SNPs written\\n")

RSCRIPT
        """
