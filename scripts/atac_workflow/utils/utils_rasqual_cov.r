# Functions
#----------------
# Covariate function for Rasqual--------------------------------

make_covariates <- function(meta_final, geno_pca, pca_points,
                            condition = c("CTL","FNF"),
                            n_pcs_max = 5,
                            output_dir = out_dir) {
  
  # Map CTL → PBS for output naming
  condition_label <- ifelse(condition == "CTL", "PBS", condition)

  # Subset metadata for the given condition
  cond_indices <- which(meta_final$Condition == condition)
  cond_samples <- meta_final[cond_indices, ]
  cond_sample_ids <- cond_samples$sampleID

  # Subset genotype PCs
  cond_geno <- geno_pca[geno_pca$sampleID %in% cond_sample_ids, ]
  cond_geno <- cond_geno[match(cond_sample_ids, cond_geno$sampleID), ]

  # Subset ATAC PCs
  cond_atac_pcs <- pca_points[cond_indices, ]

  for(i in 0:n_pcs_max) {
    message("Creating ", condition_label, " covariates with ", i, " ATAC PCs...")

    # Start with fixed covariates
    current_covars <- data.frame(
      genoPCA1 = cond_geno$genoPCA1,
      genoPCA2 = cond_geno$genoPCA2,
      genoPCA3 = cond_geno$genoPCA3,
      Sex      = cond_samples$Sex_numeric,
      Protocol = cond_samples$Protocol_batch
    )

    # Add ATAC PCs if requested
    if(i > 0) {
      atac_subset <- cond_atac_pcs[, 1:i, drop = FALSE]
      colnames(atac_subset) <- paste0("ATAC_PC", 1:i)
      current_covars <- cbind(current_covars, atac_subset)
    }

    # Save as tab-delimited text
    txt_file <- file.path(output_dir, sprintf("%s_covariates_pc%d.txt", condition_label, i))
    write.table(current_covars, txt_file,
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    # Save as binary for RASQUAL
    bin_file <- file.path(output_dir, sprintf("%s_covariates_pc%d.bin", condition_label, i))
    fbin <- file(bin_file, "wb")
    writeBin(as.double(c(as.matrix(current_covars))), fbin)
    close(fbin)
  }
}


