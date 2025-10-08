#This is MTC permutation correction
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(grid)  # for unit()
library(httpgd)

# Define parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
data_dir <- file.path(base_dir, "caQTL/data")
plot_dir <- file.path(base_dir, "caQTL/plots")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
#-------------------------------------------------
#Parameters
window_type <- "25kb"  # or 1 or 100
conditions <- c("pbs", "fnf")
pc <- paste0("pc", 0)
chromosomes <- 1:22
bf_threshold <- 0.05
#-------------------------------------------------
# permutation to combine them as one file

for (cond in conditions) {
    # Collect results from all chromosomes (eigenMT)
    chr_results <- lapply(chromosomes, function(chr) {
      file_path <- file.path(
        base_dir, "caQTL/eigenMT_perm/results",
        paste0("window_", window_type), cond, pc,
        sprintf("eigenMT_chr%s_%s_%s_perm1.txt", chr, cond, pc)
      )
      fread(file_path)
    })
    
    # Combine and filter out NULL results
    all_results <- rbindlist(chr_results[!sapply(chr_results, is.null)])
    all_results <- all_results[order(all_results$`p-value`),]
  
    output_dir <- file.path(
      base_dir, "caQTL/eigenMT_perm/combined_results",
      paste0("window_", window_type), cond, pc
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write combined results
    combined_file <- file.path(output_dir, sprintf("eigenMT_combined_%s_%s.txt", cond, pc))
    fwrite(all_results, combined_file, sep = "\t", quote = FALSE)
    
}

#---------------------------------------------------
# MTC correction using permutation results

for (cond in conditions) {
  permutation_path <- file.path(
    base_dir, "caQTL/eigenMT_perm/combined_results",
    paste0("window_", window_type), cond, pc,
    sprintf("eigenMT_combined_%s_%s.txt", cond, pc)
  )
  
  permutation <- fread(permutation_path)
  permutation <- permutation[order(permutation$`p-value`), ]
  
  ThresholdRowCounts <- floor(bf_threshold * nrow(permutation))
  thresholdvalue <- permutation$`p-value`[ThresholdRowCounts]
  
  cat("  Threshold:", sprintf("%.6e", thresholdvalue), "\n")
  
  # Apply to each chromosome
  for (chr in chromosomes) {
    
    # Load BF-adjusted RASQUAL for this chromosome
    rasqual_BFpath <- file.path(
      data_dir, "01_BF_adjusted",
      paste0("window_", window_type), cond, pc,
      sprintf("chr%s_%s_%s_%s_BFadjusted.txt", chr, cond, pc, window_type)
    )
    
    if (!file.exists(rasqual_BFpath)) {
      cat("  Chr", chr, ": file not found\n")
      next
    }
    
    rasqual <- fread(rasqual_BFpath)
    
    # Filter by threshold
    index <- which(rasqual$BF <= thresholdvalue)
    rasqual_post <- rasqual[index, ]
    
    cat("  Chr", chr, ":", nrow(rasqual_post), "significant SNPs\n")
    
    # Save per-chromosome
    output_dir <- file.path(
      data_dir, "02_MTC_final",
      paste0("window_", window_type), cond
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    output_file <- file.path(
      output_dir,
      sprintf("chr%s_%s_%s_MTCFinal.txt", chr, cond, pc)
    )
    
    fwrite(rasqual_post, output_file, sep = "\t", 
           row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}


#---- Combine all chromosome for final output---------

for (cond in conditions) {
  cat("\nProcessing:", cond, "\n")
  # Initialize combined dataframe
  Final_DF <- NULL
  
  for (chr in chromosomes) {    
    filepath <- file.path(
      data_dir, "02_MTC_final",
      paste0("window_", window_type), cond,
      sprintf("chr%s_%s_%s_MTCFinal.txt", chr, cond, pc)
    )
  
    chr_data <- fread(filepath)    
  
    if (is.null(Final_DF)) {
      Final_DF <- chr_data
    } else {
      Final_DF <- rbind(Final_DF, chr_data)
    }
    
    cat("  Chr", chr, ":", nrow(chr_data), "SNPs added\n")
  }
  
  # Save combined results
  output_dir <- file.path(
    data_dir, "02_MTC_final",
    paste0("window_", window_type), cond
  )
  
  outputfile <- file.path(
    output_dir,
    sprintf("%s_%s_AllResults_25kb_MTCFinal.txt", cond, pc)
  )
  
  fwrite(Final_DF, outputfile, sep = "\t", 
         row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  cat("\nTotal significant SNPs:", nrow(Final_DF), "\n")
  cat("Unique peaks:", length(unique(Final_DF$Feature)), "\n")
  cat("Saved to:", outputfile, "\n")
}

