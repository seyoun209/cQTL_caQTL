# ============================================================
# Analyze: caQTL + Hi-C + eQTL targets
# Find caQTLs that connect to eQTL target genes via Hi-C loops
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

# eQTL data - adjust based on your actual data structure
# Assuming pbs_eqtl and fnf_eqtl are already loaded
# Expected columns: variantID, gene_id (or gene_symbol), and ideally padj/pvalue

# ============================================================
# Helper functions (reuse from your existing code)
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
# NEW FUNCTION: Annotate caQTL-HiC-eQTL connections
# ============================================================

annotate_caqtl_hic_eqtl <- function(ld_GR, ld_anchor1_overlap, ld_anchor2_overlap,
                                    loop_gi, promoters_gr, eqtl_df,
                                    condition_name, peak_region_map) {
  
  # Prepare eQTL data
  eqtl_dt <- as.data.table(eqtl_df)
  
  # Standardize eQTL columns
  if ("gene_id" %in% names(eqtl_dt)) {
    eqtl_dt[, gene_id_clean := sub("\\..*", "", gene_id)]
  }
  
  # Get gene symbols if not present
  if (!"gene_symbol" %in% names(eqtl_dt) && "gene_id_clean" %in% names(eqtl_dt)) {
    annots <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(eqtl_dt$gene_id_clean),
      columns = c("SYMBOL"),
      keytype = "ENSEMBL"
    )
    annots <- annots %>%
      dplyr::rename(gene_id_clean = ENSEMBL, gene_symbol = SYMBOL) %>%
      distinct(gene_id_clean, .keep_all = TRUE)
    eqtl_dt <- left_join(eqtl_dt, annots, by = "gene_id_clean")
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
      
      # Check if any of these genes are eQTL targets
      for (gid_clean in gene_ids_clean) {
        # Check if this gene is an eQTL target
        eqtl_hits <- eqtl_dt[gene_id_clean == gid_clean]
        
        if (nrow(eqtl_hits) == 0) next
        
        # Take the first/best eQTL hit for this gene
        eqtl_hit <- eqtl_hits[1]
        
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
          gene_symbol    = eqtl_hit$gene_symbol,
          eqtl_variant   = if("variantID" %in% names(eqtl_hit)) eqtl_hit$variantID else NA_character_,
          eqtl_pvalue    = if("pvalue" %in% names(eqtl_hit)) eqtl_hit$pvalue else NA_real_,
          eqtl_padj      = if("padj" %in% names(eqtl_hit)) eqtl_hit$padj else NA_real_
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

cat("Annotating PBS caQTL-HiC-eQTL...\n")
pbs_caqtl_hic_eqtl <- annotate_caqtl_hic_eqtl(
  ld_GR = pbs_ld_GR,
  ld_anchor1_overlap = pbs_ld_anchor1_overlap,
  ld_anchor2_overlap = pbs_ld_anchor2_overlap,
  loop_gi = diff_loop_gi,
  promoters_gr = promoters,
  eqtl_df = pbs_eqtl,
  condition_name = "PBS",
  peak_region_map = pbs_peak_map
)

cat("Annotating FNF caQTL-HiC-eQTL...\n")
fnf_caqtl_hic_eqtl <- annotate_caqtl_hic_eqtl(
  ld_GR = fnf_ld_GR,
  ld_anchor1_overlap = fnf_ld_anchor1_overlap,
  ld_anchor2_overlap = fnf_ld_anchor2_overlap,
  loop_gi = diff_loop_gi,
  promoters_gr = promoters,
  eqtl_df = fnf_eqtl,
  condition_name = "FNF",
  peak_region_map = fnf_peak_map
)

# Combine
caqtl_hic_eqtl_all <- rbindlist(list(pbs_caqtl_hic_eqtl, fnf_caqtl_hic_eqtl), fill=TRUE)

# ============================================================
# Summary statistics
# ============================================================

cat("\n=== Summary: caQTL-HiC-eQTL connections ===\n")
print(caqtl_hic_eqtl_all[, .N, by = .(condition, region_type)])

cat("\n=== Unique genes connected ===\n")
print(caqtl_hic_eqtl_all[, .(n_genes = uniqueN(gene_id)), by = condition])

cat("\n=== Unique caQTL peaks involved ===\n")
print(caqtl_hic_eqtl_all[, .(n_peaks = uniqueN(caqtl_peak)), by = condition])

# ============================================================
# Merge with response caQTLs (optional)
# ============================================================

pbs_caqtl_hic_eqtl_response <- response_caQTL_pbs %>%
  inner_join(pbs_caqtl_hic_eqtl, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

fnf_caqtl_hic_eqtl_response <- response_caQTL_fnf %>%
  inner_join(fnf_caqtl_hic_eqtl, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

shared_caqtl_hic_eqtl <- shared_fnf %>%
  inner_join(caqtl_hic_eqtl_all, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

# ============================================================
# Save results
# ============================================================

save(pbs_caqtl_hic_eqtl, fnf_caqtl_hic_eqtl, caqtl_hic_eqtl_all,
     file = file.path(plot_data_dir, "caQTL_HiC_eQTL.RData"))

save(pbs_caqtl_hic_eqtl_response, fnf_caqtl_hic_eqtl_response, shared_caqtl_hic_eqtl,
     file = file.path(plot_data_dir, "caQTL_HiC_eQTL_response.RData"))

fwrite(caqtl_hic_eqtl_all, 
       file = file.path(plot_data_dir, "caQTL_HiC_eQTL_all.tsv"), 
       sep = "\t")


#--------------------------------------------------------------

# visualize the barplot

regionwise_summary_hic <- function(dt) {
  dt[, .(
    N = .N,
    n_genes = uniqueN(gene_id),
    n_peaks = uniqueN(caqtl_peak)
  ), by = .(condition, region_type)][order(condition, match(region_type, c("promoter","gene_body","distal")))]
}

prep_plot_df_hic <- function(sum_dt, metric = "N") {
  # metric can be "N" (connections), "n_genes", or "n_peaks"
  sum_dt[, dataset := condition]
  sum_dt[, value := get(metric)]
  sum_dt[, label := paste0(value)]
  sum_dt[]
}

plot_region_bars_hic <- function(df, title_lab, metric_label = "Number of connections") {
  color_map <- c(
    "PBS" = "#2057A7",
    "FNF" = "#F2BC40"
  )
  
  ggplot(df, aes(x = region_type, y = value, fill = dataset, group = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = label),
              position = position_dodge(width = 0.8),
              vjust = -0.3,
              size = 2, color = "black", lineheight = 0.9) +
    scale_fill_manual(values = color_map) +
    labs(x = NULL, y = metric_label, title = title_lab) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 6),
      axis.text.x = element_text(size = 4),
      axis.text.y = element_text(size = 4),
      axis.title.y = element_text(size = 6),
      axis.line.y = element_line(size = 0.2, color = "black"),
      title = element_text(size = 6),
      strip.placement = "outside",
      panel.background = element_rect(fill = "transparent", color = "transparent"),
      plot.background = element_rect(fill = "transparent", color = "transparent"),
      panel.grid = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x.bottom = element_text(size = 6),
      legend.key.size = unit(0.1, "cm"),
      legend.title = element_text(size = 4),
      legend.text = element_text(size = 4)
    )
}

# ------------------------------------------------------------
# Load caQTL-HiC-eQTL data
# ------------------------------------------------------------

load(file.path(plot_data_dir, "caQTL_HiC_eQTL.RData"))

# Summarize by region type
eqtl_hic_sum <- regionwise_summary_hic(caqtl_hic_eqtl_all)

# Prepare plotting data for different metrics
eqtl_hic_connections <- prep_plot_df_hic(eqtl_hic_sum, "N")
eqtl_hic_genes <- prep_plot_df_hic(eqtl_hic_sum, "n_genes")
eqtl_hic_peaks <- prep_plot_df_hic(eqtl_hic_sum, "n_peaks")

# Create plots
p_eqtl_hic_connections <- plot_region_bars_hic(
  eqtl_hic_connections, 
  "caQTL-HiC-eQTL connections",
  "Number of connections"
)

p_eqtl_hic_genes <- plot_region_bars_hic(
  eqtl_hic_genes, 
  "Unique eQTL genes via HiC",
  "Number of genes"
)

p_eqtl_hic_peaks <- plot_region_bars_hic(
  eqtl_hic_peaks, 
  "caQTL peaks connected to eQTL genes",
  "Number of peaks"
)

# Save plots
save(p_eqtl_hic_connections, p_eqtl_hic_genes, p_eqtl_hic_peaks,
     file = file.path(plot_data_dir, "eQTL_HiC_plots.rds"))











#--------------------------------------------------------------------------

load(file = file.path(plot_data_dir, "caQTL_HiC_eQTL.RData")) #pbs_caqtl_hic_eqtl, fnf_caqtl_hic_eqtl, caqtl_hic_eqtl_all,
load(file = file.path(plot_data_dir, "caQTL_HiC_eQTL_response.RData")) #pbs_caqtl_hic_eqtl_response, fnf_caqtl_hic_eqtl_response, shared_caqtl_hic_eqtl,


pdf(file.path(plot_dir, "highconf_eqtl_caqtl_hic_FNF.pdf"), width = 7, height = 8)
for(i in 1:nrow(fnf_caqtl_hic_eqtl_response)){
pageCreate(width = 7, height = 8, showGuides = FALSE)
plot_caqtl_multitrack_with_hic_eqtl(
  peak_id = fnf_caqtl_hic_eqtl_response$peak[i],
  snp_id = fnf_caqtl_hic_eqtl_response$snp[i],
  gene_id = fnf_caqtl_hic_eqtl_response$gene_id[i],           # ENSG ID for eQTL filtering
  gene_symbol = fnf_caqtl_hic_eqtl_response$gene_symbol[i],   # Gene symbol for display
  gene_start = fnf_caqtl_hic_eqtl_response$start[i],
  gene_end = fnf_caqtl_hic_eqtl_response$end[i],
  primary_dataset = "FNF",
  add_manhattan = TRUE,
  add_eqtl = TRUE,
  eqtl_dataset = "FNF",
  x_start = 0.5,
  y_start = 0.5,
  width = 2.5,
  height = 1.4,
  zoom_range = 400000,
  add_atac_signal = TRUE,
  add_rna_signal = TRUE
)
}
dev.off()


pdf(file.path(plot_dir, "barplot.pdf"), width = 3, height = 3)
pageCreate(width = 3, height = 3, showGuides = FALSE)
plotGG(p_eqtl_hic_connections, x = 0.5, y = 0.5, width = 2.4, height = 1)
dev.off()