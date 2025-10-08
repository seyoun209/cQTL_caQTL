library(data.table)
library(dplyr)

# Parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
chromosomes <- 1:22
conditions <- c("pbs", "fnf") 
window_kb <- 25
response_dir <- file.path(base_dir, "caQTL/data/response_qtl")
load(file.path(response_dir, "caQTL_prepdata.RData"))
load(file.path(response_dir, "rasqual_sig_caqtl.RData"))

# RASQUAL full results
rasqual_dir <- file.path(base_dir, "caQTL/rasqual_output",paste0("combined_window_", window_kb,"kb"))
pbs_rasqual <- fread(file.path(rasqual_dir, "pbs_pc0_25kb_combined.txt"))

pbs_rasqual <-  pbs_rasqual |> mutate(Condition = "pbs") |>
  dplyr::rename(peak = Feature, snp = rs_ID) |>
  mutate(Peak_SNP = paste(peak, snp, sep = "_"))

fnf_rasqual <- fread(file.path(rasqual_dir, "fnf_pc0_25kb_combined.txt"))
fnf_rasqual <-   fnf_rasqual |> mutate(Condition = "fnf") |> 
    dplyr::rename(peak = Feature, snp = rs_ID) |>
    mutate(Peak_SNP = paste(peak, snp, sep = "_"))

saveRDS(pbs_rasqual, file = file.path(response_dir, "pbs_rasqual_nominal.rds"))
saveRDS(fnf_rasqual, file = file.path(response_dir, "fnf_rasqual_nominal.rds"))

# --- Load the significant response caQTLs (FDR < 0.10) ---
sig_response_qtl <- fread(file.path(response_dir, "significant_response_caQTLs_FDR10.tsv"))


# Effectsize filter ------------------------------------------------
response_caQTL_pbs <- list()
response_caQTL_fnf <- list()

for (i in seq_len(nrow(sig_response_qtl))) {
  if (i %% 50 == 0) cat("Processing", i, "of", nrow(sig_response_qtl), "\n")

  sample <- sig_response_qtl[i, ]
  key <- sample$Peak_SNP

  # Match to PBS and FNF RASQUAL outputs
  pbs <- pbs_rasqual_sig %>% filter(Peak_SNP == key)
  fnf <- fnf_rasqual_sig %>% filter(Peak_SNP == key)
  pbs_bulk <- pbs_rasqual %>% filter(Peak_SNP == key)
  fnf_bulk <- fnf_rasqual %>% filter(Peak_SNP == key)

  # Skip if missing in both
  if (nrow(pbs) == 0 & nrow(fnf) == 0) next

  # ---- PBS condition ----
  if (nrow(pbs) > 0) {
    sample_pbs <- sample
    sample_pbs$Condition <- "pbs"
    sample_pbs$Ref <- pbs$Ref_Allele
    sample_pbs$Alt <- pbs$Alt_allele
    sample_pbs$RASQUAL_PValue_pbs <- pbs$PValue
    sample_pbs$Pi_pbs <- pbs$Effect_Size
    sample_pbs$EffectSize_pbs <- sample_pbs$Pi_pbs - 0.5

    # Compare vs FNF (if exists)
    fnf_match <- fnf_bulk %>% filter(Peak_SNP == key)
    sample_pbs$Pi_fnf <- ifelse(nrow(fnf_match) > 0, fnf_match$Effect_Size, NA)
    sample_pbs$EffectSize_fnf <- ifelse(nrow(fnf_match) > 0, fnf_match$Effect_Size - 0.5, NA)
    sample_pbs$delta_EffectSize <- sample_pbs$EffectSize_fnf - sample_pbs$EffectSize_pbs

    response_caQTL_pbs[[length(response_caQTL_pbs) + 1]] <- sample_pbs
  }

  # ---- FNF condition ----
  if (nrow(fnf) > 0) {
    sample_fnf <- sample
    sample_fnf$Condition <- "fnf"
    sample_fnf$Ref <- fnf$Ref_Allele
    sample_fnf$Alt <- fnf$Alt_allele
    sample_fnf$RASQUAL_PValue_fnf <- fnf$PValue
    sample_fnf$Pi_fnf <- fnf$Effect_Size
    sample_fnf$EffectSize_fnf <- sample_fnf$Pi_fnf - 0.5

    # Compare vs PBS (if exists)
    pbs_match <- pbs_bulk %>% filter(Peak_SNP == key)
    sample_fnf$Pi_pbs <- ifelse(nrow(pbs_match) > 0, pbs_match$Effect_Size, NA)
    sample_fnf$EffectSize_pbs <- ifelse(nrow(pbs_match) > 0, pbs_match$Effect_Size - 0.5, NA)
    sample_fnf$delta_EffectSize <- sample_fnf$EffectSize_fnf - sample_fnf$EffectSize_pbs

    response_caQTL_fnf[[length(response_caQTL_fnf) + 1]] <- sample_fnf
  }
}

# Combine into data.tables
response_caQTL_pbs <- rbindlist(response_caQTL_pbs, fill = TRUE)
response_caQTL_fnf <- rbindlist(response_caQTL_fnf, fill = TRUE)

# Save both
save(response_caQTL_pbs, file = file.path(response_dir, "response_caQTL_PBS.RData"))
save(response_caQTL_fnf, file = file.path(response_dir, "response_caQTL_FNF.RData"))


# ---- Add MAF & minor-donor estimates ----
response_caQTL_pbs <- response_caQTL_pbs %>%
  mutate(
    MAF = pmin(Allele_Frequency, 1 - Allele_Frequency),
    minor_donors = round(MAF * 2 * 21 / 2)  # ≈ per-genotype carriers
  )

response_caQTL_fnf <- response_caQTL_fnf %>%
  mutate(
    MAF = pmin(Allele_Frequency, 1 - Allele_Frequency),
    minor_donors = round(MAF * 2 * 21 / 2)
  )

# ---- High-confidence filtering ----
pbs_highconf <- response_caQTL_pbs |> 
    filter(!Peak_SNP %in% response_caQTL_fnf$Peak_SNP ) |>
    filter(
        !is.na(RASQUAL_PValue_pbs),
        RASQUAL_PValue_pbs < 0.05,
        abs(delta_EffectSize) >= 0.15,
        minor_donors >= 2) |> 
  arrange(RASQUAL_PValue_pbs)

fnf_highconf <- response_caQTL_fnf  |> 
    filter(!Peak_SNP %in% response_caQTL_pbs$Peak_SNP ) |>
    filter(
        !is.na(RASQUAL_PValue_fnf),
        RASQUAL_PValue_fnf < 0.05,
        abs(delta_EffectSize) >= 0.15,
        minor_donors >= 2
    ) |> 
    arrange(RASQUAL_PValue_fnf)


shared_fnf  <- response_caQTL_fnf |> filter(Peak_SNP %in% response_caQTL_pbs$Peak_SNP)
shared_pbs  <- response_caQTL_pbs |> filter(Peak_SNP %in% response_caQTL_fnf$Peak_SNP)

# ---- Save results ----
save(pbs_highconf, file = file.path(response_dir, "response_caQTL_PBS_highconf.RData"))
save(fnf_highconf, file = file.path(response_dir, "response_caQTL_FNF_highconf.RData"))
save(shared_pbs, shared_fnf, file = file.path(response_dir, "response_caQTL_shared.RData"))
