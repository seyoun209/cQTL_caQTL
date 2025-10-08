library(lme4)
library(lmerTest)
library(dplyr)
library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop

chrom <- as.numeric(args[1])

# Parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
chromosomes <- 1:22
conditions <- c("pbs", "fnf") 
window_kb <- 25
response_dir <- file.path(base_dir, "caQTL/data/response_qtl")
load(file.path(response_dir, "caQTL_prepdata.RData"))

# Load LD-filtered QTL pairs ------------------------------------------

#chr_file <- file.path(response_dir, paste0("chr", chrom, "_unique.tsv"))
#if (!file.exists(chr_file)) stop(paste("LD file not found:", chr_file))
#chr_qtl <- fread(chr_file)



chr_qtl <- all_caqtl_nodup %>%
  filter(grepl(paste0("^chr", chrom, "_"), peak))  


if (nrow(chr_qtl) == 0) {
 cat("No QTLs for chromosome", chrom, "\n")
 quit(save = "no")
}

# Prepare metadata -----------------------------------------------------

meta_final <- meta_final |>
  mutate(
    Condition = ifelse(Condition == "CTL", 0,
                       ifelse(Condition == "FNF", 1, NA)),
    Donor = as.factor(Donor)
  ) |>
  dplyr::select(-Sex)

# Covariate columns (exclude sampleID)
pc_cols <- setdiff(colnames(full_covar), "sampleID")
pc_formula <- paste(pc_cols, collapse = " + ")


# Function to compute minor allele count ------------------------------

calc_mac <- function(geno) {
  geno <- as.numeric(geno)
  geno <- geno[!is.na(geno)]
  if (length(geno) == 0) return(0)
  min(sum(geno), 2 * length(geno) - sum(geno))
}

# Covariate RHS for model formula -------------------------------------
pc_cols <- setdiff(colnames(full_covar), "sampleID")
pc_rhs  <- if (length(pc_cols) > 0) paste(pc_cols, collapse = " + ") else "1"


# Data frame ---------------------------------------------------------
chr_qtl <- chr_qtl |>
  mutate(
    anova_interaction_pval = as.numeric(NA),
    lrt_chisq = as.numeric(NA),
    lrt_df = as.numeric(NA),
    minor_alle_count = as.numeric(NA),
    summary_beta = as.numeric(NA),
    summary_se = as.numeric(NA),
    summary_tvalue = as.numeric(NA),
    summary_pvalue = as.numeric(NA),
    genotype_beta = as.numeric(NA),
    genotype_pval = as.numeric(NA),
    condition_beta = as.numeric(NA),
    condition_pval = as.numeric(NA)
  )

# -------------------------------
# Main LMM loop

for (i in seq_len(nrow(chr_qtl))) {

if (i %% 100 == 0) cat("Processed", i, "of", nrow(chr_qtl), "SNPs\n")

  FeatureID <- chr_qtl$peak[i]
  variantID <- chr_qtl$snp[i]
  if (!FeatureID %in% rownames(vst_counts) || !variantID %in% combined_geno_matrix$SNP) next

  # --- assemble per-feature data ---
  norm_counts_df <- vst_counts[FeatureID, ] |>
    as.numeric() |>
    (\(x) data.frame(sampleID = colnames(vst_counts), vsd = x))()

  # ---genotype data ---
  geno_row <- combined_geno_matrix |> filter(SNP == variantID)
  if (nrow(geno_row) == 0) next

  geno_samples <- setdiff(colnames(geno_row), c("CHR","SNP","(C)M","POS","COUNTED","ALT"))
  matched_samples <- intersect(colnames(vst_counts), geno_samples)
  ordered_samples <- colnames(vst_counts)[colnames(vst_counts) %in% matched_samples]
  genotype_values <- as.numeric(unlist(geno_row[, ..ordered_samples, drop = FALSE]))
  geno_df <- data.frame(sampleID = ordered_samples, genotype = genotype_values)

  # Use merge to keep only samples present in all data frames
  df <- Reduce(function(x, y) merge(x, y, by = "sampleID"),
              list(meta_final, norm_counts_df, geno_df, full_covar))

  df <- df[!is.na(df$genotype) & !is.na(df$vsd), ]

  # # --- MAC filter ---
   mac <- calc_mac(df$genotype)
   chr_qtl$minor_alle_count[i] <- mac
  # if (mac < 4) next

  # --- LMM formulas ---
  reduced_eq <- as.formula(paste("vsd ~ genotype + Condition +", pc_rhs ," + (1|Donor)"))
  full_eq    <- as.formula(paste("vsd ~ genotype + Condition + genotype:Condition + ", pc_rhs," + (1|Donor)"))

  # --- Fit models safely ---
  tryCatch({
    reduced_model <- lmer(reduced_eq, data = df, REML = FALSE)
    full_model    <- lmer(full_eq, data = df, REML = FALSE)
    res <- anova(reduced_model, full_model)

    chr_qtl$lrt_chisq[i] <- res$Chisq[2]
    chr_qtl$lrt_df[i] <- res$Df[2]
    chr_qtl$anova_interaction_pval[i] <- res$`Pr(>Chisq)`[2]

    inter_sum <- summary(full_model)$coefficients
    rownames(inter_sum) <- make.names(rownames(inter_sum))

    # --- interaction term ---
    inter_term <- grep("genotype.*Condition", rownames(inter_sum), value = TRUE)[1]
    if (!is.na(inter_term)) {
        chr_qtl$summary_beta[i]   <- inter_sum[inter_term, "Estimate"]
        chr_qtl$summary_se[i]     <- inter_sum[inter_term, "Std. Error"]
        chr_qtl$summary_tvalue[i] <- inter_sum[inter_term, "t value"]
        chr_qtl$summary_pvalue[i] <- inter_sum[inter_term, "Pr(>|t|)"]
    }

    # --- main effects ---
    if ("genotype" %in% rownames(inter_sum)) {
        chr_qtl$genotype_beta[i] <- inter_sum["genotype", "Estimate"]
        chr_qtl$genotype_pval[i] <- inter_sum["genotype", "Pr(>|t|)"]
    }
    cond_term <- grep("^Condition", rownames(inter_sum), value = TRUE)[1]
    if (!is.na(cond_term)) {
      chr_qtl$condition_beta[i] <- inter_sum[cond_term, "Estimate"]
      chr_qtl$condition_pval[i] <- inter_sum[cond_term, "Pr(>|t|)"]
    }

  }, error = function(e) {
    message("Skipping ", FeatureID, "/", variantID, " : ", e$message)
  })
}



cat("\nSummary for chromosome", chrom, ":\n")
cat("  Total QTLs tested:", nrow(chr_qtl), "\n")
#cat("  Passed MAC >= 4:", sum(!is.na(chr_qtl$minor_alle_count) & chr_qtl$minor_alle_count >= 4), "\n")
cat("  Successfully fit:", sum(!is.na(chr_qtl$anova_interaction_pval)), "\n")


out_rdata <- file.path(response_dir, paste0("chr", chrom, "_LMM_results.RData"))
out_tsv   <- gsub(".RData", ".tsv", out_rdata)
save(chr_qtl, file = out_rdata)
fwrite(chr_qtl, out_tsv, sep = "\t", quote = FALSE)