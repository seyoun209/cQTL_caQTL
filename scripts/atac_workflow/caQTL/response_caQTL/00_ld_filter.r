#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(stringr)


# ----------------------------
# Arguments and directories
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Need chromosome number as argument")

chrom <- as.numeric(args[1])
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
ld_dir <- file.path(base_dir, "caQTL/data/04_LD_results/window_25kb/ld0")
response_dir <- file.path(base_dir,"caQTL", "data", "response_qtl")
dir.create(response_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load objects
# ----------------------------
load(file.path(base_dir, "wasp", "diff_atac", "condition", "data", "00_macs2_wasp_prep_dseq2data.RData"))

# Load significant caQTLs
pbs_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results/window_25kb/pbs/pc0")
fnf_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results/window_25kb/fnf/pc0")
rasqual_dir <- file.path(base_dir, "caQTL/rasqual_output/combined_window_25kb")

pbs_caqtl <- fread(file.path(pbs_dir, "eigenMT_sigFDR_pbs_pc0.txt")) |> filter(!is.na(snp))
fnf_caqtl <- fread(file.path(fnf_dir, "eigenMT_sigFDR_fnf_pc0.txt")) |> filter(!is.na(snp))
pbs_rasqual <- fread(file.path(rasqual_dir, "pbs_pc0_25kb_combined.txt"))
fnf_rasqual <- fread(file.path(rasqual_dir, "fnf_pc0_25kb_combined.txt"))

pbs_sig <- inner_join(pbs_caqtl, pbs_rasqual, by = c("peak" = "Feature", "snp" = "rs_ID")) |> mutate(condition = "pbs")
fnf_sig <- inner_join(fnf_caqtl, fnf_rasqual, by = c("peak" = "Feature", "snp" = "rs_ID")) |> mutate(condition = "fnf")
all_caqtl <- bind_rows(pbs_sig, fnf_sig)

# ----------------------------
# LD filtering per chromosome
# ----------------------------
LD_file <- file.path(ld_dir, paste0("chr", chrom, ".LD.ld"))
LD <- fread(LD_file)
chr_caqtl <- all_caqtl |> filter(grepl(paste0("^chr", chrom, "_"), peak))
unique_peaks <- unique(chr_caqtl$peak)

all_unique_chr <- data.frame()

for (peak_id in unique_peaks) {
  peak_snps <- chr_caqtl |> filter(peak == peak_id) |> arrange(fdr, abs(Distance_From_Peak))
  ld_subset <- LD |> filter(SNP_A %in% peak_snps$snp & SNP_B %in% peak_snps$snp)
  all_unique_chr <- rbind(all_unique_chr, peak_snps[1,])

  if (nrow(peak_snps) > 1) {
    lead_snp <- peak_snps$snp[1]
    secondary <- peak_snps$snp[-1]

    ld_with_lead <- ld_subset |> 
      filter((SNP_A == lead_snp & SNP_B %in% secondary) | 
             (SNP_B == lead_snp & SNP_A %in% secondary)) |> 
      filter(R2 < 0.2)

    while (nrow(ld_with_lead) > 0) {
      snp_list <- unique(c(ld_with_lead$SNP_A, ld_with_lead$SNP_B))
      snp_list <- snp_list[snp_list != lead_snp]
      remaining <- peak_snps |> filter(snp %in% snp_list) |> arrange(fdr, abs(Distance_From_Peak))
      if (nrow(remaining) == 0) break
      all_unique_chr <- rbind(all_unique_chr, remaining[1,])
      lead_snp <- remaining$snp[1]

      if (nrow(remaining) > 1) {
        secondary <- remaining$snp[-1]
        ld_with_lead <- ld_subset |> 
          filter((SNP_A == lead_snp & SNP_B %in% secondary) | 
                 (SNP_B == lead_snp & SNP_A %in% secondary)) |> 
          filter(R2 < 0.2)
      } else break
    }
  }
}

# Remove duplicates and save
all_unique_chr$Peak_SNP <- paste(all_unique_chr$peak, all_unique_chr$snp, sep = "_")
all_unique_chr <- all_unique_chr[!duplicated(all_unique_chr$Peak_SNP), ]
fwrite(all_unique_chr, file.path(response_dir, paste0("chr", chrom, "_unique.tsv")), sep = "\t")

