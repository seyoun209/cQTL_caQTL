library(data.table)
library(dplyr)

# Parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
chromosomes <- 1:22
conditions <- c("pbs", "fnf") 
window_kb <- 25
response_dir <- file.path(base_dir, "caQTL/data/response_qtl")

#------------------------------------------------------------

all_results <- list()
for (chrom in chromosomes) {
    print(chrom)
    tsv_file <- file.path(response_dir, paste0("chr", chrom, "_LMM_results.tsv"))
  
    if (file.exists(tsv_file)) {
        chr_data <- fread(tsv_file)
        all_results[[chrom]] <- chr_data
        cat("  - Found", nrow(chr_data), "QTLs\n")
    } else {
        cat("  - File not found, skipping\n")
        }
}

# Bind all results
combined_results <- rbindlist(all_results, fill = TRUE)

# Filter out failed tests (no p-value)
combined_results_clean <- combined_results %>%
  filter(!is.na(anova_interaction_pval))

# Save combined raw results
fwrite(combined_results_clean, file.path(response_dir, "all_chromosomes_LMM_results.tsv"),
       sep = "\t", quote = FALSE)

# Multiple testing correction --------------------------------------------
combined_results_clean <- combined_results_clean %>%
  mutate(
    lrt_fdr = p.adjust(anova_interaction_pval, method = "fdr"),
    summary_fdr = p.adjust(summary_pvalue, method = "fdr")
  )

# Significant response caQTLs (FDR < 0.05)
response_caqtl_005 <- combined_results_clean %>%
  filter(lrt_fdr < 0.05) %>%
  arrange(lrt_fdr)

cat("Response caQTLs at FDR < 0.05:", nrow(response_caqtl_005), "\n")

# Significant response caQTLs (FDR < 0.1)
response_caqtl_01 <- combined_results_clean %>%
  filter(lrt_fdr < 0.1) %>%
  arrange(lrt_fdr)

cat("Response caQTLs at FDR < 0.1:", nrow(response_caqtl_01), "\n")

# Save results -----------------------------------------------------------
# All results with FDR
fwrite(combined_results_clean,
       file.path(response_dir, "all_LMM_results_with_FDR.tsv"),
       sep = "\t", quote = FALSE)

# Significant at FDR < 0.05
fwrite(response_caqtl_005,
       file.path(response_dir, "significant_response_caQTLs_FDR05.tsv"),
       sep = "\t", quote = FALSE)

# Significant at FDR < 0.1
fwrite(response_caqtl_01,
       file.path(response_dir, "significant_response_caQTLs_FDR10.tsv"),
       sep = "\t", quote = FALSE)

# Significant snp ID only
sig_response_qtl <- fread(file.path(response_dir, "significant_response_caQTLs_FDR10.tsv"))
uniqueIds <- unique(sig_response_qtl$snp)
#write.table(uniqueIds, file.path(response_dir, "sig_response_snpIDs_FDR10.txt"),
#       sep = "\t", quote = FALSE, col.names = FALSE)
uniqueIds <- data.table(V1 = unique(sig_response_qtl$snp))
uniqueIds[, chr := sub(":.*", "", V1)]
uniqueIds[, chr := gsub("^chr", "", chr)]

# Write out one file per chromosome
out_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/caQTL/data/sig_response_bychr"
dir.create(out_dir, showWarnings = FALSE)

for (c in unique(uniqueIds$chr)) {
  out_file <- file.path(out_dir, paste0("chr", c, "_sig_response_snps.txt"))
  writeLines(uniqueIds[chr == c, V1], out_file)
}