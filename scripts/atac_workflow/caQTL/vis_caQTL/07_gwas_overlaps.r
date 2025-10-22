# ============================================================
# Check GWAS overlap with caQTL-eQTL and caQTL-sQTL pairs
# ============================================================
library(data.table)
library(dplyr)
library(stringr)
library(purrr)

base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
plot_data_dir <- file.path(base_dir, "caQTL/data/plot_data")

# ------------------------------------------------------------
# Load caQTL-eQTL and caQTL-sQTL data
# ------------------------------------------------------------

load(file.path(plot_data_dir, "caQTL_eQTL_LD_annotated.Rdata"))
load(file.path(plot_data_dir, "caQTL_sQTL_LD_annotated.Rdata"))

# Filter for pairs with overlaps (eQTL_match_count > 0 or sqtl_match_count > 0)
pbs_eqtl_pairs <- pbs_eqtl_ld %>% 
  filter(eqtl_match_count > 0) %>%
  mutate(condition = "PBS", qtl_type = "eQTL")

fnf_eqtl_pairs <- fnf_eqtl_ld %>% 
  filter(eqtl_match_count > 0) %>%
  mutate(condition = "FNF", qtl_type = "eQTL")

pbs_sqtl_pairs <- pbs_sqtl_ld %>% 
  filter(sqtl_match_count > 0) %>%
  mutate(condition = "PBS", qtl_type = "sQTL")

fnf_sqtl_pairs <- fnf_sqtl_ld %>% 
  filter(sqtl_match_count > 0) %>%
  mutate(condition = "FNF", qtl_type = "sQTL")

# Combine all QTL pairs
all_qtl_pairs <- bind_rows(
  pbs_eqtl_pairs,
  fnf_eqtl_pairs,
  pbs_sqtl_pairs,
  fnf_sqtl_pairs
)

# ------------------------------------------------------------
# Parse SNP IDs
# ------------------------------------------------------------

parse_snp_id <- function(snp_id) {
  parts <- str_split(snp_id, ":", simplify = TRUE)
  data.frame(
    chrom = parts[,1],
    pos = as.numeric(parts[,2]),
    A1 = parts[,3],
    A2 = parts[,4],
    chr_pos = paste0(str_remove(parts[,1], "chr"), ":", parts[,2])
  )
}

# Parse SNPs for all QTL pairs
all_qtl_parsed <- all_qtl_pairs %>%
  bind_cols(parse_snp_id(all_qtl_pairs$snp))

cat("Total caQTL-eQTL pairs:", nrow(all_qtl_parsed[all_qtl_parsed$qtl_type == "eQTL",]), "\n")
cat("Total caQTL-sQTL pairs:", nrow(all_qtl_parsed[all_qtl_parsed$qtl_type == "sQTL",]), "\n")

# ------------------------------------------------------------
# Helper: Get LD buddies from EUR LD files
# ------------------------------------------------------------

get_eur_ld_buddies <- function(snp_id, ld_dir, r2_threshold = 0.7) {
  ld_file <- file.path(ld_dir, paste0(snp_id, ".ld"))
  if (!file.exists(ld_file)) return(data.table())
  
  ld <- tryCatch(fread(ld_file, showProgress = FALSE), error = function(e) NULL)
  if (is.null(ld) || !"R2" %in% names(ld) || !"SNP_B" %in% names(ld)) return(data.table())
  
  ld_high <- ld[R2 >= r2_threshold, .(SNP_B, R2)]
  if (nrow(ld_high) == 0) return(data.table())
  
  ld_high[, lead_snp := snp_id]
  ld_high[]
}

# ------------------------------------------------------------
# Expand QTL pairs to include EUR LD buddies
# ------------------------------------------------------------

cat("\nExpanding QTL pairs with EUR LD buddies (R² >= 0.7)...\n")
eur_ld_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05"

# Get LD buddies for all unique SNPs
unique_snps <- unique(all_qtl_parsed$snp)
cat("Processing", length(unique_snps), "unique SNPs...\n")

ld_buddies_list <- list()
pb <- txtProgressBar(min = 0, max = length(unique_snps), style = 3)
for (i in seq_along(unique_snps)) {
  ld_buddies_list[[i]] <- get_eur_ld_buddies(unique_snps[i], eur_ld_dir)
  if (i %% 100 == 0) setTxtProgressBar(pb, i)
}
close(pb)

ld_buddies_all <- rbindlist(ld_buddies_list, fill = TRUE)

if (nrow(ld_buddies_all) > 0) {
  cat("Found", nrow(ld_buddies_all), "LD buddy pairs\n")
  
  # Parse LD buddy SNP IDs
  ld_buddies_parsed <- ld_buddies_all %>%
    bind_cols(parse_snp_id(ld_buddies_all$SNP_B)) %>%
    dplyr::rename(
      buddy_snp = SNP_B,
      buddy_chrom = chrom,
      buddy_pos = pos,
      buddy_A1 = A1,
      buddy_A2 = A2
    )
  
  # Join LD buddies back to QTL pairs
  all_qtl_expanded <- all_qtl_parsed %>%
    left_join(ld_buddies_parsed, by = c("snp" = "lead_snp"), relationship = "many-to-many") %>%
    mutate(
      # Use LD buddy if available, otherwise use original SNP
      test_chrom = ifelse(!is.na(buddy_chrom), buddy_chrom, chrom),
      test_pos = ifelse(!is.na(buddy_pos), buddy_pos, pos),
      test_A1 = ifelse(!is.na(buddy_A1), buddy_A1, A1),
      test_A2 = ifelse(!is.na(buddy_A2), buddy_A2, A2),
      test_snp = ifelse(!is.na(buddy_snp), buddy_snp, snp),
      ld_R2 = ifelse(!is.na(R2), R2, 1.0)  # Original SNP has R2=1 with itself
    )
} else {
  cat("No LD buddies found, using original SNPs only\n")
  all_qtl_expanded <- all_qtl_parsed %>%
    mutate(
      test_chrom = chrom,
      test_pos = pos,
      test_A1 = A1,
      test_A2 = A2,
      test_snp = snp,
      ld_R2 = 1.0
    )
}


# ------------------------------------------------------------
# Check GWAS overlap
# ------------------------------------------------------------

# List of osteoarthritis subtypes
oa_subtypes <- c("ALLOA", "FINGER", "HAND", "HIP", "HIPKNEE", "KNEE",
                 "SPINE", "THR", "THUMB", "TJR", "TKR")

# Function to check LD overlap for each OA subtype
check_qtl_gwas_overlap <- function(oa_subtype, ancestry = "EUR", r2_threshold = 0.7) {
  
  # Read GWAS LD data
  gwas_file <- paste0("/work/users/s/e/seyoun/cQTL_caQTL/external/Hatzikotoulas_2025/hg38/", 
                      oa_subtype, "/leads/", ancestry, "_", oa_subtype, "_leads_ld_final.csv")
  
  if (!file.exists(gwas_file)) {
    warning(paste("File not found:", gwas_file))
    return(NULL)
  }
  
  gwas_ld <- fread(gwas_file)
  
  # Filter for high LD (R² >= threshold)
  gwas_ld_high <- gwas_ld %>%
    filter(ldbuddy_R2 >= r2_threshold) %>%
    # Extract chromosome and position from ldbuddy_CHR:hg38POS
    separate(`ldbuddy_CHR:hg38POS`, into = c("ld_chrom", "ld_pos"), sep = ":", remove = FALSE) %>%
    filter(!is.na(ld_pos)) %>%
    mutate(
      ld_pos = suppressWarnings(as.numeric(ld_pos)),  # Suppress NAs warning
      ld_chrom = ld_chrom
    ) %>%
    filter(!is.na(ld_pos))  # Remove rows where position couldn't be parsed
  
  # Find overlapping SNPs with explicit many-to-many relationship
  overlaps <- all_qtl_expanded %>%
    inner_join(gwas_ld_high,
               by = c("test_chrom" = "ld_chrom", "test_pos" = "ld_pos"),
               relationship = "many-to-many") %>%
    # Check allele matching
    filter(
      # Case 1: QTL alleles match GWAS LD buddy alleles in same order
      (toupper(test_A1) == toupper(ldbuddy_A1) & toupper(test_A2) == toupper(ldbuddy_A2)) |
        # Case 2: QTL alleles match GWAS LD buddy alleles in flipped order
        (toupper(test_A1) == toupper(ldbuddy_A2) & toupper(test_A2) == toupper(ldbuddy_A1))
    ) %>%
    mutate(
      oa_subtype = oa_subtype,
      ancestry = ancestry
    ) %>%
    # Keep one row per unique caQTL-GWAS pair (in case of multiple LD paths)
    group_by(snp, peak, `CHR:hg19POS`) %>%
    slice_max(order_by = ld_R2, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(snp, peak, condition, qtl_type, region_type, 
                  eqtl_annotation, eqtl_match_count, 
                  sqtl_annotation, sqtl_match_count,
                  test_snp, ld_R2,
                  oa_subtype, ancestry,
                  `CHR:hg19POS`, `CHR:hg38POS`, rsid, 
                  `Osteoarthritis phenotype`, hg19_variant,
                  metaP, ldbuddy_R2, ldbuddy_A1, ldbuddy_A2)
  
  return(overlaps)
}

# Check overlaps for all OA subtypes
cat("\nChecking GWAS overlaps...\n")
all_qtl_gwas_overlaps <- map_dfr(oa_subtypes, ~check_qtl_gwas_overlap(.x, ancestry = "EUR"))

# ------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("GWAS OVERLAP SUMMARY\n")
cat(rep("=", 60), "\n\n", sep = "")

# Overall summary
overall_summary <- all_qtl_gwas_overlaps %>%
  group_by(qtl_type) %>%
  summarise(
    n_qtl_pairs = n_distinct(paste(snp, peak)),
    n_qtl_snps = n_distinct(snp),
    n_gwas_leads = n_distinct(`CHR:hg19POS`),
    n_rsid = n_distinct(rsid),
    n_oa_subtypes = n_distinct(oa_subtype),
    .groups = "drop"
  )

cat("Overall summary by QTL type:\n")
print(overall_summary)

# By condition and QTL type
condition_summary <- all_qtl_gwas_overlaps %>%
  group_by(condition, qtl_type) %>%
  summarise(
    n_qtl_pairs = n_distinct(paste(snp, peak)),
    n_qtl_snps = n_distinct(snp),
    n_gwas_leads = n_distinct(`CHR:hg19POS`),
    n_rsid = n_distinct(rsid),
    .groups = "drop"
  ) %>%
  arrange(qtl_type, condition)

cat("\n\nSummary by condition and QTL type:\n")
print(condition_summary)

# By OA subtype
oa_summary <- all_qtl_gwas_overlaps %>%
  group_by(oa_subtype, qtl_type) %>%
  summarise(
    n_qtl_pairs = n_distinct(paste(snp, peak)),
    n_qtl_snps = n_distinct(snp),
    n_gwas_leads = n_distinct(`CHR:hg19POS`),
    n_rsid = n_distinct(rsid),
    .groups = "drop"
  ) %>%
  arrange(qtl_type, oa_subtype)

cat("\n\nSummary by OA subtype and QTL type:\n")
print(oa_summary)

# By region type
region_summary <- all_qtl_gwas_overlaps %>%
  group_by(region_type, qtl_type) %>%
  summarise(
    n_qtl_pairs = n_distinct(paste(snp, peak)),
    n_qtl_snps = n_distinct(snp),
    n_gwas_leads = n_distinct(`CHR:hg19POS`),
    .groups = "drop"
  ) %>%
  arrange(qtl_type, region_type)

cat("\n\nSummary by region type and QTL type:\n")
print(region_summary)

# ------------------------------------------------------------
# Check novelty
# ------------------------------------------------------------

coloc_table <- fread("/work/users/s/e/seyoun/cQTL_caQTL/external/Hatzikotoulas_2025/gwas_sup19.csv")

# Add novelty check
all_qtl_gwas_with_novelty <- all_qtl_gwas_overlaps %>%
  mutate(
    # Check if this GWAS variant appears in colocalization table
    known_coloc_lead = hg19_variant %in% coloc_table$lead.SNP,
    known_gwas_signal = hg19_variant %in% coloc_table$indep.gwas.snp,
    
    # Overall novelty status
    novelty_status = case_when(
      known_coloc_lead ~ "Known colocalization lead",
      known_gwas_signal ~ "Known GWAS signal (tested in coloc)",
      TRUE ~ "Potentially novel"
    )
  )

# Novelty summary
novelty_summary <- all_qtl_gwas_with_novelty %>%
  group_by(qtl_type, novelty_status) %>%
  summarise(
    n_qtl_pairs = n_distinct(paste(snp, peak)),
    n_qtl_snps = n_distinct(snp),
    n_gwas_leads = n_distinct(`CHR:hg19POS`),
    .groups = "drop"
  )

cat("\n\nNovelty summary:\n")
print(novelty_summary)

# Potentially novel findings
novel_findings <- all_qtl_gwas_with_novelty %>%
  filter(novelty_status == "Potentially novel")

cat("\n\nPotentially novel findings:\n")
cat("  eQTL:", nrow(novel_findings[novel_findings$qtl_type == "eQTL",]), "overlaps\n")
cat("  sQTL:", nrow(novel_findings[novel_findings$qtl_type == "sQTL",]), "overlaps\n")


save(all_qtl_gwas_overlaps, all_qtl_gwas_with_novelty, 
     overall_summary, condition_summary, oa_summary, region_summary, novelty_summary,
     file = file.path(plot_data_dir, "caQTL_QTL_GWAS_overlaps.RData"))

fwrite(all_qtl_gwas_overlaps, 
       file = file.path(plot_data_dir, "caQTL_QTL_GWAS_overlaps.tsv"),
       sep = "\t")

fwrite(all_qtl_gwas_with_novelty, 
       file = file.path(plot_data_dir, "caQTL_QTL_GWAS_overlaps_with_novelty.tsv"),
       sep = "\t")





table_phenotype_variant <- all_qtl_gwas_overlaps %>%
  mutate(peak_snp = paste(peak, snp, sep = "-")) %>%
  group_by(rsid, condition, qtl_type) %>%
  summarise(
    n_phenotypes = n_distinct(`Osteoarthritis phenotype`),
    phenotypes = paste(unique(`Osteoarthritis phenotype`), collapse = ","),
    peak_snp_pairs = paste(unique(peak_snp), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_phenotypes >= 1) |> arrange(qtl_type) |> as.data.frame()
