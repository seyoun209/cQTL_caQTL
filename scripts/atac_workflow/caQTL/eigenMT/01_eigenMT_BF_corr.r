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
window_type <- "25kb"  # or 1, 10,25, or 100
conditions <- c("pbs", "fnf")
pcs <- paste0("pc", 0:10)
chromosomes <- 1:22
fdr_threshold <- 0.05
bonf_threshold <- 0.05
#-------------------------------------------------

for (cond in conditions) {
  for (pc in pcs) {
    cat("Processing:", cond, pc, "\n")  
    # Process each chromosome separately
    lapply(chromosomes, function(chr) {
      
      # Load eigenMT (lead SNP + TESTS)
      eigenmt_path <- file.path(
        base_dir, "caQTL/eigenMT/results",
        paste0("window_", window_type), cond, pc,
        sprintf("eigenMT_chr%s_%s_%s.txt", chr, cond, pc)
      )
      
      eigenmt <- fread(eigenmt_path)
      colnames(eigenmt)[2] <- "Feature"
      
      # Load RASQUAL (all SNPs)
      rasqual_path <- file.path(
        base_dir, "caQTL/rasqual_output",
        paste0("filtered_window_", window_type), pc,
        sprintf("%s_chr%s_%s_%s_filtered.txt", cond, chr, pc, window_type)
      )
      
      rasqual <- fread(rasqual_path)
      
      # Check match
      rasq_feats <- unique(rasqual$Feature)
      eigen_feats <- unique(eigenmt$Feature)
      
      if (length(rasq_feats) != length(eigen_feats)) {
        cat(sprintf("    Mismatch: RASQUAL=%d, eigenMT=%d peaks\n",
                    length(rasq_feats), length(eigen_feats)))
      }
      
      # Merge - only take TESTS from eigenMT
      rasqual <- rasqual %>%
        left_join(eigenmt %>% select(Feature, TESTS), by = "Feature")
      
      # Recalculate BF for ALL SNPs
      rasqual$PValue <- as.numeric(rasqual$PValue)
      rasqual$TESTS <- as.numeric(rasqual$TESTS)
      rasqual$BF <- rasqual$PValue * rasqual$TESTS
      rasqual$BF <- ifelse(rasqual$BF > 1, 1, rasqual$BF)
      
      # Save per-chromosome BF-adjusted results
      output_dir <- file.path(
        data_dir, "01_BF_adjusted",
        paste0("window_", window_type), cond, pc
      )
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      output_file <- file.path(output_dir,sprintf("chr%s_%s_%s_%s_BFadjusted.txt", chr, cond, pc, window_type))
      fwrite(rasqual, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
      cat(sprintf("    Saved %d SNPs\n", nrow(rasqual)))
      
      return(NULL)
    })
  }
}



