# ============================================================
# Analyze: caQTL + Hi-C + sQTL targets
# Find caQTLs that connect to sQTL target genes via Hi-C loops
# ============================================================

library(data.table)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(InteractionSet)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ------------------------------------------------------------
# Load processed objects
# ------------------------------------------------------------
base_dir       <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
response_dir   <- file.path(base_dir, "caQTL/data/response_qtl")
plot_data_dir  <- file.path(base_dir, "caQTL/data/plot_data")
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(response_dir, "rasqual_sig_caqtl.RData"))        # pbs_rasqual_sig, fnf_rasqual_sig
load(file.path(response_dir, "response_caQTL_PBS.RData"))       # response_caQTL_pbs
load(file.path(response_dir, "response_caQTL_FNF.RData"))       # response_caQTL_fnf
load(file.path(response_dir, "response_caQTL_shared.RData"))    # shared_pbs, shared_fnf

# Hi-C data
hic_processed_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/caQTL_processing/CQTL_AI/output/hic/00_data"
load(file.path(hic_processed_dir, "diff_loop_gi.Rdata"))        # diff_loop_gi

# Gene annotation
txdb <- loadDb("/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb")
txdb_genes <- genes(txdb)
promoters  <- promoters(txdb_genes, upstream = 2500, downstream = 500)
exons_gr   <- exons(txdb)
genes_gr   <- genes(txdb)

# sQTL data - assuming pbs_sQTL and fnf_sQTL are already loaded
# Expected columns: var_id, phe_id, SYMBOL, and ideally padj/pvalue

# ============================================================
# Helper functions (reuse from existing code)
# ============================================================

make_caQTL_GRanges <- function(peak_vec) {
  peaks <- unique(na.omit(peak_vec))
  split <- str_split_fixed(peaks, "_", 3)
  GRanges(seqnames = split[,1],
          ranges   = IRanges(start = as.numeric(split[,2]), end = as.numeric(split[,3])))
}

get_caqtl_ld_variants <- function(caPeaks, ld_dir, r2_threshold = 0.7) {
  result_list <- list()
  for (i in seq_len(nrow(caPeaks))) {
    snp_id <- caPeaks$snp[i]
    ld_file <- file.path(ld_dir, paste0(snp_id, ".ld"))
    if (!file.exists(ld_file)) next
    ld <- tryCatch(fread(ld_file), error = function(e) NULL)
    if (is.null(ld)) next
    hi <- ld[R2 > r2_threshold]
    if (nrow(hi) > 0) {
      hi$lead_snp <- snp_id
      hi$peak <- caPeaks$peak[i]
      result_list[[length(result_list)+1]] <- hi
    }
  }
  if (length(result_list)==0) return(data.table())
  rbindlist(result_list, fill = TRUE)
}

collapse_ld_variants <- function(ld_dt) {
  if (nrow(ld_dt)==0) return(ld_dt)
  ld_dt %>%
    group_by(peak, lead_snp) %>%
    arrange(desc(R2)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    as.data.table()
}

convert_ld_to_granges <- function(ld_dt) {
  if (nrow(ld_dt)==0) return(GRanges())
  snp_split <- str_split_fixed(ld_dt$SNP_B, ":", 4)
  GRanges(seqnames = snp_split[,1],
          ranges = IRanges(start = as.numeric(snp_split[,2]), width = 1),
          ld_snp = ld_dt$SNP_B,
          lead_caqtl = ld_dt$lead_snp,
          peak = ld_dt$peak,
          R2 = ld_dt$R2)
}

label_peak_region_type <- function(peaks_char, promoters_gr, genes_gr) {
  gr <- make_caQTL_GRanges(peaks_char)
  prom_ov <- countOverlaps(gr, promoters_gr) > 0
  gene_ov <- countOverlaps(gr, genes_gr) > 0
  region <- ifelse(prom_ov, "promoter",
                   ifelse(gene_ov, "gene_body", "distal"))
  data.table(peak = peaks_char, region_type = region)
}

# ============================================================
# Classify caQTL peaks by region type
# ============================================================

pbs_peak_map <- label_peak_region_type(pbs_rasqual_sig$peak, promoters, genes_gr)
fnf_peak_map <- label_peak_region_type(fnf_rasqual_sig$peak, promoters, genes_gr)

# ============================================================
# Get LD variants and convert to GRanges
# ============================================================

pbs_ld <- get_caqtl_ld_variants(pbs_rasqual_sig, 
                                ld_dir="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05",
                                r2_threshold=0.7)
fnf_ld <- get_caqtl_ld_variants(fnf_rasqual_sig,
                                ld_dir="/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05",
                                r2_threshold=0.7)

pbs_ld_collapsed <- collapse_ld_variants(pbs_ld)
fnf_ld_collapsed <- collapse_ld_variants(fnf_ld)

pbs_ld_GR <- convert_ld_to_granges(pbs_ld_collapsed)
fnf_ld_GR <- convert_ld_to_granges(fnf_ld_collapsed)

# Get Hi-C anchors
anchor1 <- anchors(diff_loop_gi, type="first")
anchor2 <- anchors(diff_loop_gi, type="second")

# Find overlaps between LD variants and Hi-C anchors
pbs_ld_anchor1_overlap <- findOverlaps(pbs_ld_GR, anchor1)
pbs_ld_anchor2_overlap <- findOverlaps(pbs_ld_GR, anchor2)
fnf_ld_anchor1_overlap <- findOverlaps(fnf_ld_GR, anchor1)
fnf_ld_anchor2_overlap <- findOverlaps(fnf_ld_GR, anchor2)

# ============================================================
# NEW FUNCTION: Annotate caQTL-HiC-sQTL connections
# ============================================================

annotate_caqtl_hic_sqtl <- function(ld_GR, ld_anchor1_overlap, ld_anchor2_overlap,
                                    loop_gi, promoters_gr, sqtl_df,
                                    condition_name, peak_region_map) {
  
  # Prepare sQTL data
  sqtl_dt <- as.data.table(sqtl_df)
  
  # Standardize sQTL columns - extract gene_id from phe_id if needed
  if ("phe_id" %in% names(sqtl_dt) && !"gene_id" %in% names(sqtl_dt)) {
    # phe_id format is typically like "ENSG00000123456:E001:E002"
    sqtl_dt[, gene_id := sub(":.*", "", phe_id)]
  }
  
  if ("gene_id" %in% names(sqtl_dt)) {
    sqtl_dt[, gene_id_clean := sub("\\..*", "", gene_id)]
  }
  
  # Get gene symbols if not present
  if (!"SYMBOL" %in% names(sqtl_dt) && "gene_id_clean" %in% names(sqtl_dt)) {
    annots <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(sqtl_dt$gene_id_clean),
      columns = c("SYMBOL"),
      keytype = "ENSEMBL"
    )
    annots <- annots %>%
      dplyr::rename(gene_id_clean = ENSEMBL, SYMBOL = SYMBOL) %>%
      distinct(gene_id_clean, .keep_all = TRUE)
    sqtl_dt <- left_join(sqtl_dt, annots, by = "gene_id_clean")
  }
  
  results <- list()
  
  # Internal helper function
  add_hits <- function(hits_obj, ld_anchor_label, gene_anchor_label) {
    if (length(hits_obj) == 0) return(invisible(NULL))
    qh <- queryHits(hits_obj)
    sh <- subjectHits(hits_obj)
    
    for (i in seq_along(qh)) {
      ld_idx   <- qh[i]
      loop_idx <- sh[i]
      peak_name <- ld_GR$peak[ld_idx]
      
      # Get region type
      region_type <- peak_region_map[peak == peak_name, region_type]
      if (length(region_type) == 0) region_type <- NA
      
      lead_snp <- ld_GR$lead_caqtl[ld_idx]
      ld_snp   <- ld_GR$ld_snp[ld_idx]
      r2       <- ld_GR$R2[ld_idx]
      
      # Get opposite anchor (gene side)
      gene_anchor <- anchors(loop_gi, type = ifelse(grepl("1$", ld_anchor_label),
                                                    "second", "first"))[loop_idx]
      
      # Find promoters at gene anchor
      promoter_ov <- findOverlaps(gene_anchor, promoters_gr)
      if (length(promoter_ov) == 0) next
      gene_ids <- promoters_gr$gene_id[subjectHits(promoter_ov)]
      gene_ids_clean <- sub("\\..*", "", gene_ids)
      
      # Check if any of these genes have sQTLs
      for (gid_clean in gene_ids_clean) {
        # Check if this gene has sQTLs
        sqtl_hits <- sqtl_dt[gene_id_clean == gid_clean]
        
        if (nrow(sqtl_hits) == 0) next
        
        # Take the first/best sQTL hit for this gene
        sqtl_hit <- sqtl_hits[1]
        
        results[[length(results) + 1]] <<- data.table(
          condition      = condition_name,
          caqtl_peak     = peak_name,
          region_type    = region_type,
          lead_snp       = lead_snp,
          ld_variant     = ld_snp,
          R2             = r2,
          loop_index     = loop_idx,
          ld_at_anchor   = ld_anchor_label,
          gene_at_anchor = gene_anchor_label,
          gene_id        = gid_clean,
          gene_symbol    = sqtl_hit$SYMBOL,
          sqtl_variant   = if("var_id" %in% names(sqtl_hit)) sqtl_hit$var_id else NA_character_,
          sqtl_phe_id    = if("phe_id" %in% names(sqtl_hit)) sqtl_hit$phe_id else NA_character_,
          sqtl_pvalue    = if("pvalue" %in% names(sqtl_hit)) sqtl_hit$pvalue else NA_real_,
          sqtl_padj      = if("padj" %in% names(sqtl_hit)) sqtl_hit$padj else NA_real_
        )
      }
    }
  }
  
  # Process both anchor overlaps
  add_hits(ld_anchor1_overlap, "anchor1", "anchor2")
  add_hits(ld_anchor2_overlap, "anchor2", "anchor1")
  
  if (length(results) == 0) return(data.table())
  out <- rbindlist(results)
  
  # Deduplicate
  setorder(out, -R2)
  out <- out[!duplicated(out[, .(caqtl_peak, gene_id, loop_index)]), ]
  
  return(out[])
}

# ============================================================
# Run for PBS and FNF
# ============================================================

cat("Annotating PBS caQTL-HiC-sQTL...\n")
pbs_caqtl_hic_sqtl <- annotate_caqtl_hic_sqtl(
  ld_GR = pbs_ld_GR,
  ld_anchor1_overlap = pbs_ld_anchor1_overlap,
  ld_anchor2_overlap = pbs_ld_anchor2_overlap,
  loop_gi = diff_loop_gi,
  promoters_gr = promoters,
  sqtl_df = pbs_sQTL,
  condition_name = "PBS",
  peak_region_map = pbs_peak_map
)

cat("Annotating FNF caQTL-HiC-sQTL...\n")
fnf_caqtl_hic_sqtl <- annotate_caqtl_hic_sqtl(
  ld_GR = fnf_ld_GR,
  ld_anchor1_overlap = fnf_ld_anchor1_overlap,
  ld_anchor2_overlap = fnf_ld_anchor2_overlap,
  loop_gi = diff_loop_gi,
  promoters_gr = promoters,
  sqtl_df = fnf_sQTL,
  condition_name = "FNF",
  peak_region_map = fnf_peak_map
)

# Nothing Found. 

# Combine
caqtl_hic_sqtl_all <- rbindlist(list(pbs_caqtl_hic_sqtl, fnf_caqtl_hic_sqtl), fill=TRUE)

# ============================================================
# Summary statistics
# ============================================================

cat("\n=== Summary: caQTL-HiC-sQTL connections ===\n")
print(caqtl_hic_sqtl_all[, .N, by = .(condition, region_type)])

cat("\n=== Unique genes connected ===\n")
print(caqtl_hic_sqtl_all[, .(n_genes = uniqueN(gene_id)), by = condition])

cat("\n=== Unique caQTL peaks involved ===\n")
print(caqtl_hic_sqtl_all[, .(n_peaks = uniqueN(caqtl_peak)), by = condition])

# ============================================================
# Merge with response caQTLs (optional)
# ============================================================

pbs_caqtl_hic_sqtl_response <- response_caQTL_pbs %>%
  inner_join(pbs_caqtl_hic_sqtl, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

fnf_caqtl_hic_sqtl_response <- response_caQTL_fnf %>%
  inner_join(fnf_caqtl_hic_sqtl, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

shared_caqtl_hic_sqtl <- shared_fnf %>%
  inner_join(caqtl_hic_sqtl_all, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

# ============================================================
# Save results
# ============================================================

save(pbs_caqtl_hic_sqtl, fnf_caqtl_hic_sqtl, caqtl_hic_sqtl_all,
     file = file.path(plot_data_dir, "caQTL_HiC_sQTL.RData"))

save(pbs_caqtl_hic_sqtl_response, fnf_caqtl_hic_sqtl_response, shared_caqtl_hic_sqtl,
     file = file.path(plot_data_dir, "caQTL_HiC_sQTL_response.RData"))

fwrite(caqtl_hic_sqtl_all, 
       file = file.path(plot_data_dir, "caQTL_HiC_sQTL_all.tsv"), 
       sep = "\t")

