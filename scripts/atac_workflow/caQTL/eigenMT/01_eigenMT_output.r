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

# Main processing loop
for (cond in conditions) {
  for (pc in pcs) {
    
    cat("Processing:", cond, pc, "\n")
    
    # Collect results from all chromosomes
    chr_results <- lapply(chromosomes, function(chr) {
      file_path <- file.path(
        base_dir, "caQTL/eigenMT/results",
        paste0("window_", window_type), cond, pc,
        sprintf("eigenMT_chr%s_%s_%s.txt", chr, cond, pc)
      )
      
      if (!file.exists(file_path)) {
        cat("  Warning: Missing", basename(file_path), "\n")
        return(NULL)
      }
      
      fread(file_path)
    })
    
    # Combine and filter out NULL results
    all_results <- rbindlist(chr_results[!sapply(chr_results, is.null)])
    
    if (nrow(all_results) == 0) {
      cat("  No results found, skipping\n\n")
      next
    }
    
    # Calculate FDR from Bonferroni-corrected p-values
    all_results <- all_results %>%
      arrange("p-value") %>%
      mutate(fdr = p.adjust(BF, method = "fdr"))
    
    # Identify significant associations
    sig_results <- all_results %>% filter(fdr < fdr_threshold)
    n_sig <- nrow(sig_results)
    
    # Get maximum nominal p-value among significant hits
    pval_cutoff <- if (n_sig > 0) max(sig_results$"p-value") else NA_real_
    
    # Create output directory
    output_dir <- file.path(
      base_dir, "caQTL/eigenMT/combined_results",
      paste0("window_", window_type), cond, pc
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write combined results
    combined_file <- file.path(output_dir, sprintf("eigenMT_combined_%s_%s.txt", cond, pc))
    fwrite(all_results, combined_file, sep = "\t", quote = FALSE)
    
    # Write p-value cutoff
    cutoff_file <- file.path(output_dir, sprintf("%s_%s_MaxPvalueCutoff_%02d.txt", 
                                                  cond, pc, fdr_threshold * 100))
    fwrite(data.table(max_pvalue = pval_cutoff), cutoff_file, sep = "\t")
    
    cat(sprintf("  Total: %d | Significant (FDR < %.2f): %d\n\n", 
                nrow(all_results), fdr_threshold, n_sig))
  }
}

cat("All eigenMT results processed successfully.\n")
#----------------------------------------------------------------
# Visualization of significant SNP counts across PCs
# Count significant associations for each condition/PC combination
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
    n_significant = {
      if (file.exists(combined_file)) {
        results <- fread(combined_file, select = "fdr")
        sum(results$fdr < fdr_threshold, na.rm = TRUE)
      } else {
        NA_integer_
      }
    }
  ) %>%
  ungroup() %>%
  select(condition, pc, n_significant)

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

# Find maximum counts per condition for annotation
max_summary <- summary_df %>%
  group_by(type) %>%
  slice_max(n_significant, n = 1, with_ties = TRUE) %>%
  mutate(
    color = case_when(
      type == "PBS" ~ "#2057A7",
      type == "FN-f" ~ "#F2BC40",
      TRUE ~ "#999999"
    )
  ) %>%
  ungroup()

# Save summary
summary_output <- file.path(
  base_dir, "caQTL/eigenMT/combined_results",
  paste0("window_", window_type),
  sprintf("summary_significant_caQTLs_FDR%02d.txt", fdr_threshold * 100)
)
fwrite(summary_df, summary_output, sep = "\t")


counts_dotPlot <- ggplot(summary_df, aes(x = pc, y = n_significant, color = type)) +
  geom_point(size = 2) +
  labs(x = "PC", y = "Number of significant caQTLs", color = "Type") +
  scale_color_manual(values = c("PBS" = "#2057A7", "FN-f" = "#F2BC40")) +
  geom_text(data = max_summary, aes(x = pc, y = n_significant + 100,
                                   label = paste0(pc, ": ", n_significant)),
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
save(counts_dotPlot, file = file.path(plot_dir, "qc", "counts_dotplot.RData"))
ggsave(filename = file.path(plot_dir, "qc", "counts_dotplot.pdf"),
       plot = counts_dotPlot, width = 7, height = 4, units = "in")
print(counts_dotPlot)
