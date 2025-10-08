library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

#-------------------------------------------------
# Input directories
#-------------------------------------------------
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
combined_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", "window_25kb")
geno_dir <- file.path(base_dir, "geno")  # has pbs_geno and fnf_geno dirs
plot_dir <- file.path(base_dir, "caQTL/plots/top_hits")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
conditions <- c("pbs", "fnf")
pcs <- paste0("pc", 0)
fdr_threshold <- 0.05
n_top <- 20  # number of hits to visualize per condition

#-------------------------------------------------
# Helper function: read eigenMT + filter
#-------------------------------------------------
get_top_hits <- function(cond, pc, n_top=20) {
  file_path <- file.path(combined_dir, cond, pc, sprintf("eigenMT_combined_%s_%s.txt", cond, pc))
  if (!file.exists(file_path)) return(NULL)
  df <- fread(file_path)
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>% mutate(fdr = p.adjust(`p-value`, method="fdr"))
  df <- df %>% filter(fdr < fdr_threshold)
  if (nrow(df) == 0) return(NULL)
  
  df %>% arrange(`p-value`) %>% head(n_top)
}

#-------------------------------------------------
# Collect top hits for PBS and FN-f
#-------------------------------------------------
top_hits <- list()

for (cond in conditions) {
  cond_hits <- lapply(pcs, function(pc) get_top_hits(cond, pc, n_top)) %>%
    bind_rows() %>% mutate(condition = cond)
  if (!is.null(cond_hits)) top_hits[[cond]] <- cond_hits
}

all_hits <- bind_rows(top_hits)
fwrite(all_hits, file.path(plot_dir, "top_hits_summary.txt"), sep="\t")

#-------------------------------------------------
# Function to load genotype from .traw
#-------------------------------------------------
load_genotypes <- function(chrom, cond) {
  traw_file <- file.path(geno_dir, paste0(cond, "_geno"), sprintf("chr%s.%s_geno.traw", chrom, cond))
  if (!file.exists(traw_file)) return(NULL)
  geno <- fread(traw_file)
  # Format: CHR SNP ID ... then sample genotypes
  geno_long <- geno %>%
    pivot_longer(cols=starts_with("CQTL"), names_to="sampleID", values_to="genotype")
  return(geno_long)
}

#-------------------------------------------------
# Example: plot top hits (PBS + FN-f combined)
#-------------------------------------------------
for (i in 1:nrow(all_hits)) {
  hit <- all_hits[i,]
  chrom <- hit$chr
  var_id <- hit$var_id
  phe_id <- hit$phe_id
  cond <- hit$condition
  
  # Load genotype
  geno_df <- load_genotypes(chrom, cond)
  if (is.null(geno_df)) next
  geno_df <- geno_df %>% filter(SNP == var_id)
  
  # Load expression (replace with your own PSI/counts matrix)
  # For demo: assume `expr_matrix` exists with rownames=peaks (phe_id), cols=sampleID
  expr_row <- expr_matrix[phe_id,]
  expr_df <- data.frame(sampleID=names(expr_row), expr=as.numeric(expr_row))
  
  plot_df <- left_join(geno_df, expr_df, by="sampleID") %>%
    mutate(genotype = factor(genotype, levels=c(0,1,2)))
  
  p <- ggplot(plot_df, aes(x=genotype, y=expr, fill=genotype)) +
    geom_boxplot() +
    geom_jitter(width=0.2, alpha=0.5) +
    labs(title=paste(cond, phe_id, var_id, sep=" | "),
         x="Genotype", y="Expression (normalized)") +
    theme_minimal()
  
  ggsave(file.path(plot_dir, sprintf("%s_chr%s_%s_%s.pdf", cond, chrom, phe_id, var_id)),
         plot=p, width=4, height=4)
}
