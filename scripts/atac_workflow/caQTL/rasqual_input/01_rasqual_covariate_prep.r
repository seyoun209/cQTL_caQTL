# Rasqual - Covariate Preparation
# Seyoun Byun
# Date: 2025-09-24
# RASQUAL requires covariate in the bin format (WASP)   
#----------------------------------------------------------------
# Load Libraries
library(data.table)
library(dplyr)
library(tibble)
library(base)
library(DESeq2)

#----------------------------------------------------------------
# Load Data
source("/work/users/s/e/seyoun/cQTL_caQTL/scripts/atac_workflow/utils/utils_diff_atac.r")
source("/work/users/s/e/seyoun/cQTL_caQTL/scripts/atac_workflow/utils/utils_rasqual_cov.r")
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
geno_dir <- file.path(base_dir,"geno")
out_dir <-file.path(base_dir, "caQTL","rasqual_input","covar")
n_pcs_max <- 10

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


#----------------------------------------------------------------
# 1. atac-seq PC calculation
load(file.path(base_dir,"wasp","diff_atac","condition","data","01_macs2_wasp_cqn_normFactors.RData")) #cqn.counts, cqnNormFactors,GCcontent,peakwidths
load(file.path(base_dir,"wasp","diff_atac","condition","data","00_macs2_wasp_prep_dseq2data.RData")) # se, macs2_gr, counts_macs2, meta_final

atac_dds_rasqual <- DESeqDataSetFromMatrix(countData = counts_macs2,
                                  colData = meta_final,
                                  rowRanges=macs2_gr,
                                  design = ~Donor)

normalizationFactors(atac_dds_rasqual) <- cqnNormFactors
keep <- rowSums(counts(atac_dds_rasqual) >= 10)  >= ceiling(nrow(colData(se))*0.10)
atac_dds_rasqual <- atac_dds_rasqual[keep,]

# Get VST transformed data
vsd_donor <- vst(atac_dds_rasqual, blind=TRUE)
vst_counts = assay(vsd_donor)

# Calculate PCA - this step was missing!
pca_rasqual <- prcomp(t(vst_counts))

# Correlation plot
# Create PCA data frame
pca_df <- as.data.frame(pca_rasqual$x)
pca_vsd_df <- merge(pca_df, colData(vsd_donor), by = "row.names", all.x = TRUE)
ctl_plot <- corr_plot(pca_vsd_df, condition = "CTL", title = "CTL Correlations")
fnf_plot <- corr_plot(pca_vsd_df, condition = "FNF", title = "FNF Correlations")

# Calculate variance explained
var_explained <- pca_rasqual$sdev^2 / sum(pca_rasqual$sdev^2)
message("Variance explained by first 42 PCs: ",
        paste(round(var_explained[1:42]*100, 2), collapse="%, "), "%")
pca_points <- as.data.frame(pca_rasqual$x[,1:42])

#----------------------------------------------------------------
# 2.clean up genotype PCs
## Read PBS genotype PCA
pbs_geno_pca <- fread(file.path(geno_dir,"pbs_geno","02_pca","cqtl.eigenvec"))
pbs_geno_pca <- pbs_geno_pca %>%
  dplyr::select(-V2) %>%
  dplyr::rename(sampleID = V1) %>%
  rename_with(~ paste0("genoPCA", 1:20), .cols = V3:V22)

## Read FNF genotype PCA
fnf_geno_pca <- fread(file.path(geno_dir,"fnf_geno","02_pca","cqtl.eigenvec"))
fnf_geno_pca <- fnf_geno_pca %>%
  dplyr::select(-V2) %>%
  dplyr::rename(sampleID = V1) %>%
  rename_with(~ paste0("genoPCA", 1:20), .cols = V3:V22)


#----------------------------------------------------------------
# 3. Clean up metadata
## Add numeric versions of key covariates
meta_final <- meta_final %>%
  mutate(
    Sex_numeric = as.numeric(Sex == "F"),  # 0=Male, 1=Female
    Protocol_batch = as.numeric(factor(ATACLibrarySubmissionDate))  # Protocol batch as numeric factor
  )

#---------------------------------------------------------------
# 4. Generate covariates
#PBS
make_covariates(meta_final, pbs_geno_pca, pca_points,
                condition = "CTL", n_pcs_max = 10,
                output_dir = out_dir)

# FNF
make_covariates(meta_final, fnf_geno_pca, pca_points,
                condition = "FNF", n_pcs_max = 10,
                output_dir = out_dir)


#-----------------------------------------------------------------------

# 5. Save mapping information for reference
# Make sure indices and sample IDs exist
pbs_indices <- which(meta_final$Condition == "CTL")
fnf_indices <- which(meta_final$Condition == "FNF")

pbs_sample_ids <- meta_final$sampleID[pbs_indices]
fnf_sample_ids <- meta_final$sampleID[fnf_indices]

# Create mapping with PBS instead of CTL
sample_mapping <- data.frame(
  condition = c(rep("PBS", length(pbs_indices)),
                rep("FNF", length(fnf_indices))),
  sample_id = c(pbs_sample_ids, fnf_sample_ids),
  sample_index = c(pbs_indices, fnf_indices)
)

# Save
write.table(sample_mapping,
            file.path(out_dir, "sample_mapping.txt"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t"
            )


# PBS mapping
pbs_mapping <- data.frame(
  condition    = rep("PBS", length(pbs_indices)),
  sample_id    = pbs_sample_ids,
  sample_index = seq_len(length(pbs_indices)) 
)

# FNF mapping
fnf_mapping <- data.frame(
  condition    = rep("FNF", length(fnf_indices)),
  sample_id    = fnf_sample_ids,
  sample_index = seq_len(length(fnf_indices))
)

# Save separately
write.table(pbs_mapping,
            file.path(out_dir, "sample_mapping_PBS.txt"),
            row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")

write.table(fnf_mapping,
            file.path(out_dir, "sample_mapping_FNF.txt"),
            row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")



message("All covariates created successfully. Files saved to: ", out_dir)
message("File naming: pbs_covariates_pc0.bin through pbs_covariates_pc10.bin")
message("File naming: fnf_covariates_pc0.bin through fnf_covariates_pc10.bin")