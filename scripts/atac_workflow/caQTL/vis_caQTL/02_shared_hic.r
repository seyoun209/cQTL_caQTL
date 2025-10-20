library(data.table)
library(dplyr)
library(GenomeInfoDb)
library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggvenn)
library(scales)
library(colorspace)
library(ggtext)
library(ggforce)
library(httpgd)
library(plotgardener)
library(colorspace)
library(mariner)
library(InteractionSet)

#------------------------------------------------------------
# Load data
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
response_dir <- file.path(base_dir, "caQTL/data/response_qtl")
plot_dir <- file.path(base_dir, "caQTL/plots")
plot_data_dir <- file.path(base_dir, "caQTL/data/plot_data")
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)



load(file.path(response_dir, "caQTL_prepdata.RData")) #vst_counts, peak_info,  all_rasqual_qtl_duplicated,all_caqtl_nodup,  all_caqtl, combined_geno_matrix, full_covar, meta_final

load(file.path(response_dir, "response_caQTL_PBS.RData")) #response_caQTL_pbs
load(file.path(response_dir, "response_caQTL_FNF.RData")) # response_caQTL_fnf
load(file.path(response_dir, "response_caQTL_PBS_highconf.RData"))
load(file.path(response_dir, "response_caQTL_FNF_highconf.RData"))
load(file.path(response_dir, "response_caQTL_shared.RData")) #shared_pbs, shared_fnf
pbs_rasqual <- readRDS(file.path(response_dir, "pbs_rasqual_nominal.rds"))
fnf_rasqual <- readRDS(file.path(response_dir, "fnf_rasqual_nominal.rds"))

load(file.path(response_dir, "rasqual_sig_caqtl.RData")) # pbs_rasqual_sig, fnf_rasqual_sig

# Load eigenMT results
chromosomes <- 1:22
conditions <- c("pbs", "fnf")
window_kb <- 25

# eigenMT significant directories
pbs_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "pbs", "pc0")
fnf_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "fnf", "pc0")


# eQTL data
eqtl_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/eqtl/"
pbs_eqtl <- fread(paste0(eqtl_dir,"CTL_PEER_k20_genoPC_perm1Mb_sig_rsID_signalRanges.csv"))
fnf_eqtl <- fread(paste0(eqtl_dir,"FNF_PEER_k22_genoPC_perm1Mb_sig_rsID_signalRanges.csv"))

# RNA gene expression data
de_genes <- fread("/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/quant/de_genes_results.csv")

# sQTL data
response_pbs_results <- readRDS("/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")
# Subset significant sQTL only
pbs_sQTL <- response_pbs_results %>%
  dplyr::filter(rank  == 0)

fnf_sQTL <- response_fnf_results %>%
  dplyr::filter(rank  == 0)

# hic data
hic_processed_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/caQTL_processing/CQTL_AI/output/hic/00_data"
hic_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/hic_caQTL/dietJuicer/output"
load(file.path(hic_processed_dir, "diff_loop_gi.Rdata")) # diff_loop_gi, this is gi interaction clean
load(file.path(hic_processed_dir, "all_sig_loops.Rdata")) # sig_loops, this is sig loops padj < 0.1 and no log2 foldchange
load(file.path(hic_processed_dir, "fnf_gained_loops.Rdata")) # fnf_gained only sig for the fnf
load(file.path(hic_processed_dir, "pbs_gained_loops.Rdata")) # pbs_gained, fnf_lost for the fnf but I labeld wrong as the pbs gained(not gained but fnf lost

#make TXDB for the v49
#This is only one time
#txdb <- makeTxDbFromGFF(file="/proj/phanstiel_lab/Reference/human/hg38//annotations/gencode.v49.primary_assembly.annotation.gtf")
#saveDb(x=txdb, file = "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb")
txdb <- loadDb("/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb")

txdb_genes <- genes(txdb)
promoters <- promoters(txdb_genes, upstream =2500, downstream =500)

#---- sig SNPS 
pbs_caSNP <- unique(pbs_rasqual_sig$snp)
fnf_caSNP <- unique(fnf_rasqual_sig$snp)
total_caQTL_unique <- unique(c(pbs_caSNP, fnf_caSNP))

write.table(total_caQTL_unique,
    file.path(base_dir, "geno/sigcaQTLtotal.list"),
            col.names = F,
            row.names=F,
            quote=F)

uniqueIds <- data.table(V1 = unique(total_caQTL_unique))
uniqueIds[, chr := sub(":.*", "", V1)]
uniqueIds[, chr := gsub("^chr", "", chr)]

# Write out one file per chromosome
out_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/caQTL/data/sig_response_bychr"
#dir.create(out_dir, showWarnings = FALSE)

for (c in unique(uniqueIds$chr)) {
  out_file <- file.path(out_dir, paste0("chr", c, "_sig_response_snps.txt"))
  writeLines(uniqueIds[chr == c, V1], out_file)
}

#----------------------------------------------------------------------------------------
# Overlapping hic-loop anchors from the loops -----------------------------
# Analyze1: caQTL + Hic loop only 
#HIC loop overlaps------------------------------------------------------------------------
# I have exteended the caqtl peak by +/- 25kb
pbs_peak_GR <- make_caQTL_GRanges_expanded(pbs_rasqual_sig$peak, window_kb = window_kb)
fnf_peak_GR <- make_caQTL_GRanges_expanded(fnf_rasqual_sig$peak, window_kb = window_kb)
anchor1 <- anchors(diff_loop_gi, type = "first")
anchor2 <- anchors(diff_loop_gi, type = "second")

# PBS overlaps
pbs_peak_anchor1_overlap <- findOverlaps(pbs_peak_GR, anchor1)
pbs_peak_anchor2_overlap <- findOverlaps(pbs_peak_GR, anchor2)

# FNF overlaps
fnf_peak_anchor1_overlap <- findOverlaps(fnf_peak_GR, anchor1)
fnf_peak_anchor2_overlap <- findOverlaps(fnf_peak_GR, anchor2)


# Getting unique caQTL or hic loop hits
pbs_caqtl_hits <- unique(c(queryHits(pbs_peak_anchor1_overlap), queryHits(pbs_peak_anchor2_overlap)))
fnf_caqtl_hits <- unique(c(queryHits(fnf_peak_anchor1_overlap), queryHits(fnf_peak_anchor2_overlap)))

pbs_loops_hits <- unique(c(subjectHits(pbs_peak_anchor1_overlap),subjectHits(pbs_peak_anchor2_overlap)))
fnf_loops_hits <- unique(c(subjectHits(fnf_peak_anchor1_overlap),subjectHits(fnf_peak_anchor2_overlap)))

pbs_caqtl_overlaps <- pbs_rasqual_sig[pbs_caqtl_hits, ]
fnf_caqtl_overlaps <- fnf_rasqual_sig[fnf_caqtl_hits, ]
pbs_loops_overlaps <- diff_loop_gi[pbs_loops_hits, ]
fnf_loops_overlaps <- diff_loop_gi[fnf_loops_hits, ]

# Run LD lookup R2 > 0.7 ( I need to change it to own sample)
pbs_ld_variants <- get_caqtl_ld_variants(pbs_caqtl_overlaps, 
                                          ld_dir = "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05",
                                          r2_threshold = 0.7)

fnf_ld_variants <- get_caqtl_ld_variants(fnf_caqtl_overlaps, 
                                          ld_dir = "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05",
                                          r2_threshold = 0.7)

save(pbs_ld_variants, fnf_ld_variants, file = file.path(plot_data_dir, "caQTL_hic_overlaps_LD_variants.RData"))




# Check LD > 0.7 snps overlps with anchors of hic 

pbs_ld_GR <- convert_ld_to_granges(pbs_ld_variants)
fnf_ld_GR <- convert_ld_to_granges(fnf_ld_variants)

pbs_ld_anchor1_overlap <- findOverlaps(pbs_ld_GR, anchor1)
pbs_ld_anchor2_overlap <- findOverlaps(pbs_ld_GR, anchor2)

fnf_ld_anchor1_overlap <- findOverlaps(fnf_ld_GR, anchor1)
fnf_ld_anchor2_overlap <- findOverlaps(fnf_ld_GR, anchor2)

#----------------------------------------------------------

# --- 1) classify peaks as promoter vs distal --------------------------------

pbs_peak_map <- label_peak_region_type(pbs_rasqual_sig$peak, promoters)
fnf_peak_map <- label_peak_region_type(fnf_rasqual_sig$peak, promoters)

# --- 2) collapse LD buddies per peak ----------------------------------------

pbs_ld_variants_collapsed <- collapse_ld_variants(pbs_ld_variants)
fnf_ld_variants_collapsed <- collapse_ld_variants(fnf_ld_variants)

pbs_ld_GR <- convert_ld_to_granges(pbs_ld_variants_collapsed)
fnf_ld_GR <- convert_ld_to_granges(fnf_ld_variants_collapsed)

# --- 3) overlaps of (collapsed) LD with anchors ------------------------------

anchor1 <- anchors(diff_loop_gi, type = "first")
anchor2 <- anchors(diff_loop_gi, type = "second")

pbs_ld_anchor1_overlap <- findOverlaps(pbs_ld_GR, anchor1)
pbs_ld_anchor2_overlap <- findOverlaps(pbs_ld_GR, anchor2)

fnf_ld_anchor1_overlap <- findOverlaps(fnf_ld_GR, anchor1)
fnf_ld_anchor2_overlap <- findOverlaps(fnf_ld_GR, anchor2)

# --- 4) annotate LONG-RANGE (distal) caQTL → loop → promoter → DE gene ------

# ensure de_genes has 'gene_id', 'log2FoldChange', 'padj'
de_genes_clean <- as.data.table(de_genes) |> mutate(gene_id_noVer = gsub("\\..*", "", gene_id))
cat("Annotating PBS long-range...\n")
pbs_longrange <- annotate_longrange_ld_loop_genes(
  ld_GR                = pbs_ld_GR,
  ld_anchor1_overlap   = pbs_ld_anchor1_overlap,
  ld_anchor2_overlap   = pbs_ld_anchor2_overlap,
  loop_gi              = diff_loop_gi,
  promoters_gr         = promoters,
  de_genes_df          = de_genes_clean,
  condition_name       = "PBS",
  peak_region_map      = pbs_peak_map,
  padj_thresh          = 0.05
)

cat("Annotating FNF long-range...\n")
fnf_longrange <- annotate_longrange_ld_loop_genes(
  ld_GR                = fnf_ld_GR,
  ld_anchor1_overlap   = fnf_ld_anchor1_overlap,
  ld_anchor2_overlap   = fnf_ld_anchor2_overlap,
  loop_gi              = diff_loop_gi,
  promoters_gr         = promoters,
  de_genes_df          = de_genes_clean,
  condition_name       = "FNF",
  peak_region_map      = fnf_peak_map,
  padj_thresh          = 0.05
)
save(pbs_longrange, fnf_longrange, file = file.path(plot_data_dir, "longrange_caqtl_loop_DE.RData"))
# --- 5) tidy outputs + quick summaries --------------------------------------

longrange_all <- rbindlist(list(pbs_longrange, fnf_longrange), fill = TRUE)

# optional: add a clean Ensembl ID without version and (if you want) a symbol map
longrange_all[, gene_id_clean := sub("\\..*", "", gene_id)]

# summary counts
summary_tbl <- longrange_all[, .(
  n_links = .N,
  n_unique_peaks = uniqueN(caqtl_peak),
  n_unique_genes = uniqueN(gene_id_clean),
  median_abs_log2fc = median(abs(de_log2fc), na.rm = TRUE)
), by = condition][order(condition)]

print(summary_tbl)


# check whether it's also response caQTL

pbs_longrange_merged <- response_caQTL_pbs %>%
  inner_join(pbs_longrange, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

# For FNF
fnf_longrange_merged <- response_caQTL_fnf %>%
  inner_join(fnf_longrange, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

#shared

shared_longrange_merged <- shared_fnf %>%
  inner_join(longrange_all, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))


# --- 6) save ----------------------------------------------------------------

saveRDS(longrange_all, file.path(plot_data_dir, "longrange_caqtl_loop_DE_collapsed.rds"))
fwrite(longrange_all, file.path(plot_data_dir, "longrange_caqtl_loop_DE_collapsed.tsv"), sep = "\t")









# #----------------------------------------------------------------------------------------
# # Analyze 2: caQTL + eQTL-----------------------------------------------------------------

# #----------------------------------------------------------------------------------------
# # Analyze 3: caQTL + sQTL-----------------------------------------------------------------


# #----------------------------------------------------------------------------------------
# # Analyze 4: caQTL + Hic +eQTL------------------------------------------------------------

# #----------------------------------------------------------------------------------------
# # Analyze 5: caQTL + Hic + sQTL-----------------------------------------------------------