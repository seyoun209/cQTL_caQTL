# This is for the linear mixed model caQTL  preparation
# Libraries
library(DESeq2)
library(data.table)
library(dplyr)
library(stringr)
#------------------------------------------------------------
# Parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
chromosomes <- 1:22
conditions <- c("pbs", "fnf")  # Your conditions
window_kb <- 25

# eigenMT significant directories
pbs_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "pbs", "pc0")
fnf_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "fnf", "pc0")

# RASQUAL full results
rasqual_dir <- file.path(base_dir, "caQTL/rasqual_output",paste0("combined_window_", window_kb,"kb"))

# LD directory
ld_dir <- file.path(base_dir, "caQTL/data/04_LD_results/window_25kb/ld0")

# Geno directory
pbs_geno_dir <- file.path(base_dir, "geno/pbs_geno/04_recode/")
fnf_geno_dir <- file.path(base_dir, "geno/fnf_geno/04_recode/")

# meta file  (se, macs2_gr, counts_macs2, meta_final)
load(file.path(base_dir,"wasp","diff_atac","condition","data","00_macs2_wasp_prep_dseq2data.RData"))

response_qtl_dir <- file.path(base_dir, "caQTL/data/response_qtl")

#------------------------------------------------------------
# eigenMT significant QTLs------------------------
# PBS
pbs_caqtl <- fread(file.path(pbs_dir, "eigenMT_sigFDR_pbs_pc0.txt"))
pbs_snps <- pbs_caqtl %>% filter(!is.na(snp))
pbs_caqtl$condition <- "pbs"

# FNf
fnf_caqtl <- fread(file.path(fnf_dir, "eigenMT_sigFDR_fnf_pc0.txt"))
fnf_snps <- fnf_caqtl %>% filter(!is.na(snp))
fnf_caqtl$condition <- "fnf"

all_caqtl <- rbind(pbs_caqtl, fnf_caqtl)


# RASQUAL FULL results---------------------------
pbs_rasqual <- fread(file.path(rasqual_dir, "pbs_pc0_25kb_combined.txt"))
fnf_rasqual <- fread(file.path(rasqual_dir, "fnf_pc0_25kb_combined.txt"))

pbs_rasqual_sig <- inner_join(pbs_caqtl, pbs_rasqual, by = c("peak" = "Feature", "snp" = "rs_ID")) |>
  mutate(Condition = "pbs") |> 
  mutate(Peak_SNP = paste(peak, snp, sep = "_"))
fnf_rasqual_sig <- inner_join(fnf_caqtl, fnf_rasqual, by = c("peak" = "Feature", "snp" = "rs_ID")) |>
mutate(Condition = "fnf") |> 
mutate(Peak_SNP = paste(peak, snp, sep = "_"))

save(pbs_rasqual_sig, fnf_rasqual_sig, file = file.path(response_qtl_dir, "rasqual_sig_caqtl.RData"))

# Combine all significant caQTLs
all_caqtl <- rbind(pbs_rasqual_sig, fnf_rasqual_sig)
#all_caqtl$Peak_SNP <- paste(all_caqtl$peak, all_caqtl$snp, sep = "_")


#Making a duplicated SNPs object, can use to reference caQTLs significant in BOTH conditions
all_rasqual_qtl_duplicated <- all_caqtl[duplicated(all_caqtl[, "Peak_SNP"]),] # n= 9456
#Unique Peak-SNP pairings list significant in EITHER condition; those with duplicates noted above
all_caqtl_nodup <- all_caqtl[!duplicated(all_caqtl[, "Peak_SNP"]),] # n= 21871


# ATAC Counts-------------------------------------------------

# Read the saved count tables
load(file.path(base_dir,"caQTL","data","response_qtl","vst_counts.RData")) # vst_counts also peak_info

# GENOTYPES-------------------------------------------------

read_geno_matrix <- function(geno_dir, chromosomes, snp_list, prefix) {
  files <- paste0(geno_dir, "/chr", chromosomes, ".", prefix, "_geno.traw")
  geno_list <- lapply(files, fread)
  geno_matrix <- bind_rows(geno_list)
  geno_matrix %>% filter(SNP %in% snp_list)
}
pbs_geno_matrix_final <- read_geno_matrix(pbs_geno_dir, chromosomes, all_caqtl$snp, "pbs")
fnf_geno_matrix_final <- read_geno_matrix(fnf_geno_dir, chromosomes, all_caqtl$snp, "fnf")

# Define key columns for merging
key_cols <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")

# Merge PBS and FNF genotype matrices by key columns
combined_geno_matrix <- full_join(
  pbs_geno_matrix_final, fnf_geno_matrix_final,
  by = key_cols
)

# covariates------------------------------------------------
pbs_cov <- fread(file.path(base_dir, "caQTL/rasqual_input/covar/pbs_covariates_pc0.txt"))
fnf_cov <- fread(file.path(base_dir, "caQTL/rasqual_input/covar/fnf_covariates_pc0.txt"))

# Sample info ------------------------------------------------
meta_final <- meta_final %>%
  mutate(
    Sex_numeric = as.numeric(Sex == "F"),  # 0=Male, 1=Female
    Protocol_batch = as.numeric(factor(ATACLibrarySubmissionDate))  # Numeric batch
  )

# Get sample IDs for each condition
ctl_samples <- meta_final %>% filter(Condition == "CTL") %>% pull(sampleID)
fnf_samples <- meta_final %>% filter(Condition == "FNF") %>% pull(sampleID)

# Combine covariates with sample IDs
covar_pbs_final <- bind_cols(pbs_cov, sampleID = ctl_samples)
covar_fnf_final <- bind_cols(fnf_cov, sampleID = fnf_samples)

# Final covariate table
full_covar <- bind_rows(covar_pbs_final, covar_fnf_final)


# LD filtering with R2 < 0.2-----------------------------------------
# all_caqtl_unique <- data.frame()

# for (chrom in chromosomes) {
#     print(chrom)
#     LD_file <- paste0(ld_dir, "/chr", chrom, ".LD.ld")
#     LD <- fread(LD_file)

# # Filter caQTLs for this chromosome
#     chr_caqtl <- all_caqtl %>% filter(grepl(paste0("^chr", chrom, "_"), peak))

# # Get unique peaks
#     unique_peaks <- unique(chr_caqtl$peak)
  
# # LD filtering per peak
#     for (peak_id in unique_peaks) {
#         peak_snps <- chr_caqtl |> filter(peak == peak_id) |> 
#             arrange(fdr, abs(Distance_From_Peak))

# # Filter LD to only SNPs for this peak
#         ld_subset <- LD |> filter(SNP_A %in% peak_snps$snp & SNP_B %in% peak_snps$snp)

# # Add first lead SNP (best FDR)
#         all_caqtl_unique <- rbind(all_caqtl_unique, peak_snps[1,])
    
#         if (nrow(peak_snps) > 1) {
# # Lead SNP
#         lead_snp <- peak_snps$snp[1]
#         secondary_snps <- peak_snps$snp[2:nrow(peak_snps)]
      
# # Find SNPs NOT in LD with lead (R2 < 0.2)
#         ld_with_lead <- ld_subset %>%
#         filter((SNP_A == lead_snp & SNP_B %in% secondary_snps) | 
#                (SNP_B == lead_snp & SNP_A %in% secondary_snps)) %>%
#         filter(R2 < 0.2) 
      
#         # Iteratively add independent SNPs
#         while (nrow(ld_with_lead) > 0) {
#             snp_list <- unique(c(ld_with_lead$SNP_A, ld_with_lead$SNP_B))
#             snp_list <- snp_list[snp_list != lead_snp]

#             remaining <- peak_snps |> filter(snp %in% snp_list) |>
#                 arrange(fdr, abs(Distance_From_Peak))
        
#         if (nrow(remaining) == 0) break
        
#         # Add next independent SNP
#         all_caqtl_unique <- rbind(all_caqtl_unique, remaining[1,])
        
#         # Update lead
#         lead_snp <- remaining$snp[1]
        
#         if (nrow(remaining) > 1) {
#           secondary_snps <- remaining$snp[2:nrow(remaining)]
          
#           ld_with_lead <- ld_subset %>%
#             filter((SNP_A == lead_snp & SNP_B %in% secondary_snps) | 
#                    (SNP_B == lead_snp & SNP_A %in% secondary_snps)) %>%
#             filter(R2 < 0.2)
#         } else {
#           break
#         }
#       }
#     }
#   }
# }

# all_caqtl_unique$Peak_SNP <- paste(all_caqtl_unique$peak, all_caqtl_unique$snp, sep = "_")
# all_caqtl_unique <- all_caqtl_unique[!duplicated(all_caqtl_unique$Peak_SNP), ]


save(vst_counts, peak_info,  all_rasqual_qtl_duplicated,all_caqtl_nodup,  all_caqtl,
     combined_geno_matrix, full_covar, meta_final,
     file = file.path(response_qtl_dir, "caQTL_prepdata.RData")) #all_caqtl_unique

response_dir <- file.path(base_dir, "caQTL/data/response_qtl")
files <- list.files(response_dir, pattern = "chr[0-9]+_unique.tsv", full.names = TRUE)
all_caqtl_unique <- rbindlist(lapply(files, fread))
all_caqtl_unique <- all_caqtl_unique[!duplicated(all_caqtl_unique$Peak_SNP), ]
save(all_caqtl_unique, file = file.path(response_dir, "all_caqtl_unique.RData"))