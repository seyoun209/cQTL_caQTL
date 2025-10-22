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
pcs <- paste0("pc", 0:10)
chromosomes <- 1:22
fdr_threshold <- 0.05
bonf_threshold <- 0.05
#-------------------------------------------------
for (cond in conditions) {
  for (pc_cov in pcs) {
    cat("Processing:", cond, pc_cov, "\n")

    allchr_results <- NULL

    # Process each chromosome
    for (chr in chromosomes) {
      file_path <- file.path(
        base_dir, "caQTL/eigenMT/results",
        paste0("window_", window_type), cond, pc_cov,
        sprintf("eigenMT_chr%d_%s_%s.txt", chr, cond, pc_cov)
      )

      if (!file.exists(file_path)) {
        cat("  Warning: File not found -", file_path, "\n")
        next
      }

      eigen_results <- fread(file_path)
      allchr_results <- rbind(allchr_results, eigen_results)
    }

    # Skip if no data
    if (is.null(allchr_results) || nrow(allchr_results) == 0) {
      cat("  No results found for condition:", cond, "and PC:", pc_cov, "\n")
      next
    }

    # Ensure numeric columns
    allchr_results$`p-value` <- as.numeric(allchr_results$`p-value`)
    allchr_results$BF <- as.numeric(allchr_results$BF)

    # Sort by nominal p-value
    allchr_results <- allchr_results %>% arrange(`p-value`)

    # Compute global FDR using BF column
    allchr_results$fdr <- p.adjust(allchr_results$BF, method = "fdr", n = nrow(allchr_results))

    # Significant SNPs (FDR < 0.05)
    sig_snps <- allchr_results %>% filter(fdr < 0.05)
    pval_max <- if (nrow(sig_snps) > 0) max(sig_snps$`p-value`) else NA

    # Output directory
    output_dir <- file.path(
      base_dir, "caQTL/eigenMT/combined_results",
      paste0("window_", window_type), cond, pc_cov
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Write outputs
    fwrite(allchr_results, file.path(output_dir,
           sprintf("eigenMT_combined_%s_%s.txt", cond, pc_cov)),
           sep = "\t", quote = FALSE)

    fwrite(sig_snps, file.path(output_dir,
           sprintf("eigenMT_sigFDR_%s_%s.txt", cond, pc_cov)),
           sep = "\t", quote = FALSE)

    write.table(pval_max, file.path(output_dir,
                sprintf("%s_%s_MaxPvalueCutoff_05.txt", cond, pc_cov)),
                col.names = TRUE, row.names = FALSE, quote = FALSE)

    # Log
    cat(sprintf("  Total results: %d | FDR < 0.05: %d | Max P = %s\n\n",
                nrow(allchr_results), nrow(sig_snps),
                ifelse(is.na(pval_max), "NA", format(pval_max, digits = 3))))
  }
}



# # Main processing loop
# for (cond in conditions) {
#   for (pc in pcs) {
#     cat("Processing:", cond, pc, "\n")
    
#     chr_results <- lapply(chromosomes, function(chr) {
#       file_path <- file.path(
#         base_dir, "caQTL/eigenMT/results",
#         paste0("window_", window_type), cond, pc,
#         sprintf("eigenMT_chr%s_%s_%s.txt", chr, cond, pc)
#       )
#       if (file.exists(file_path)) fread(file_path) else NULL
#     })
    
#     all_results <- rbindlist(chr_results[!sapply(chr_results, is.null)])
#     if (nrow(all_results) == 0) {
#       cat("  No results found, skipping\n\n")
#       next
#     }
    
#     # Bonferroni
#     n_tests <- nrow(all_results)
#     genome_bonf_threshold <- bonf_threshold / n_tests
#     sig_bonf <- all_results %>% filter(BF < genome_bonf_threshold)
#     n_sig_bonf <- nrow(sig_bonf)
    
#     # FDR
#     all_results <- all_results %>% mutate(fdr = p.adjust(BF, method = "fdr"))
#     sig_fdr <- all_results %>% filter(fdr < fdr_threshold)
#     n_sig_fdr <- nrow(sig_fdr)
    
#     # Log results
#     cat(sprintf("  Bonferroni-significant: %d (%.2f%%)\n", n_sig_bonf, 100*n_sig_bonf/n_tests))
#     cat(sprintf("  FDR-significant: %d (%.2f%%)\n", n_sig_fdr, 100*n_sig_fdr/n_tests))
    
#     # Write outputs
#     output_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_type), cond, pc)
#     dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
#     fwrite(all_results, file.path(output_dir, sprintf("eigenMT_combined_%s_%s.txt", cond, pc)), sep="\t", quote=FALSE)
#     fwrite(sig_bonf, file.path(output_dir, sprintf("eigenMT_sigBonf_%s_%s.txt", cond, pc)), sep="\t", quote=FALSE)
#     fwrite(sig_fdr, file.path(output_dir, sprintf("eigenMT_sigFDR_%s_%s.txt", cond, pc)), sep="\t", quote=FALSE)
    
#     summary_table <- data.table(
#       total_tests = n_tests,
#       bonf_cutoff = genome_bonf_threshold,
#       n_sig_bonf = n_sig_bonf,
#       n_sig_fdr = n_sig_fdr
#     )
#     fwrite(summary_table, file.path(output_dir, sprintf("summary_%s_%s.txt", cond, pc)), sep="\t")
#   }
# }

#----------------------------------------------------------------
# Visualization of significant SNP counts across PCs
summary_df <- expand.grid(
  condition = conditions,
  pc = pcs,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    combined_file = file.path(
      base_dir, "caQTL/eigenMT/combined_results",
      paste0("window_", window_type), condition, pc,
      sprintf("eigenMT_combined_%s_%s.txt", condition, pc)
    ),
    n_bonferroni = {
      if (file.exists(combined_file)) {
        results <- fread(combined_file)
        n_tests <- nrow(results)
        genome_bonf <- 0.05 / n_tests
        sum(results$BF < genome_bonf, na.rm = TRUE)
      } else {
        NA_integer_
      }
    },
    n_fdr = {
      if (file.exists(combined_file)) {
        results <- fread(combined_file, select = "fdr")
        sum(results$fdr < 0.05, na.rm = TRUE)
      } else {
        NA_integer_
      }
    }
  ) %>%
  ungroup() %>%
  select(condition, pc, n_bonferroni, n_fdr)

# Format for plotting
summary_df <- summary_df %>%
  mutate(
    type = toupper(condition),
    pc = factor(pc, levels = pcs),
    type = recode(type, "FNF" = "FN-f")
  ) %>%
  mutate(
    type = factor(type, levels = c("PBS", "FN-f"))
  )

# Save summary
summary_output <- file.path(
  base_dir, "caQTL/eigenMT/combined_results",
  paste0("window_", window_type),
  "summary_both_Bonferroni_FDR.txt"
)
fwrite(summary_df, summary_output, sep = "\t")


#-----------------------------------------------------------------
summary_df  <- fread(summary_output)

# Find maximum counts per condition for annotation
max_summary <- summary_df %>%
  group_by(type) %>%
  slice_max(n_fdr, n = 1, with_ties = TRUE) %>%
  mutate(
    color = case_when(
      type == "PBS" ~ "#2057A7",
      type == "FN-f" ~ "#F2BC40",
      TRUE ~ "#999999"
    )
  ) %>%
  ungroup()


counts_dotPlot <- ggplot(summary_df, aes(x = pc, y = n_fdr, color = type)) +
  geom_point(size = 2) +
  labs(x = "PC", y = "Number of significant caQTLs", color = "Type") +
  scale_color_manual(values = c("PBS" = "#2057A7", "FN-f" = "#F2BC40")) +
  geom_text(data = max_summary, aes(x = pc, y = n_fdr + 100,
                                   label = paste0(pc, ": ", n_fdr)),
            color = max_summary$color, vjust = 0, size = 3.5) +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.position.inside = c(0.9, 0.9),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8),
        axis.ticks.x = element_line(color = "black", linewidth = 0.1),
        axis.text.x = element_text(color = "black", size = 8)
  )

# Save the plot to a file and display it
dir.create(file.path(plot_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
saveRDS(counts_dotPlot, file = file.path(plot_dir, "qc", paste0(window_type, "_counts_dotplot.rds")))
ggsave(filename = file.path(plot_dir, "qc", paste0(window_type, "_counts_dotplot.pdf")),
       plot = counts_dotPlot, width = 7, height = 4, units = "in")
print(counts_dotPlot)



# Make a table for the  1kb, 10kb, 25kb, 100kb