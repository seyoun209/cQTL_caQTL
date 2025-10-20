# ============================================================
# Analyze 1: caQTL + Hi-C + DE genes (using relaxed promoters)
# ============================================================

library(data.table)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(InteractionSet)
library(mariner)
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
load(file.path(response_dir, "caQTL_prepdata.RData"))           # meta, peak_info, etc.

# Hi-C data
hic_processed_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/caQTL_processing/CQTL_AI/output/hic/00_data"
load(file.path(hic_processed_dir, "diff_loop_gi.Rdata"))        # diff_loop_gi

# Gene annotation
txdb <- loadDb("/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb")
txdb_genes <- genes(txdb)
promoters  <- promoters(txdb_genes, upstream = 2500, downstream = 500)   # ← relaxed promoter window
exons_gr   <- exons(txdb)
genes_gr   <- genes(txdb)

# DE gene table
de_genes <- fread("/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/quant/de_genes_results.csv")
de_genes_clean <- as.data.table(de_genes)[, gene_id_noVer := gsub("\\..*", "", gene_id)]

annots <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = de_genes_clean$gene_id_noVer,
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)

annots <- annots %>%
  dplyr::rename(gene_id_noVer = ENSEMBL) %>%
  distinct(gene_id_noVer, .keep_all = TRUE)

de_genes_clean <- left_join(de_genes_clean, annots, by = "gene_id_noVer")

# ============================================================
# helper functions
# ============================================================

make_caQTL_GRanges <- function(peak_vec) {
  peaks <- unique(na.omit(peak_vec))
  split <- str_split_fixed(peaks, "_", 3)
  GRanges(seqnames = split[,1],
          ranges   = IRanges(start = as.numeric(split[,2]), end = as.numeric(split[,3])))
}

make_caQTL_GRanges_expanded <- function(peak_vec, window_kb = 25) {
  peaks <- unique(na.omit(peak_vec))
  split <- str_split_fixed(peaks, "_", 3)
  chr   <- split[,1]
  start <- as.numeric(split[,2])
  end   <- as.numeric(split[,3])
  center <- floor((start + end)/2)
  w <- window_kb*1000
  GRanges(seqnames = chr,
          ranges = IRanges(start = pmax(center - w, 1), end = center + w))
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
# classify caQTL peaks by region type
# ============================================================

pbs_peak_map <- label_peak_region_type(pbs_rasqual_sig$peak, promoters, genes_gr)
fnf_peak_map <- label_peak_region_type(fnf_rasqual_sig$peak, promoters, genes_gr)

# ============================================================
# LD, overlap, and annotation
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

anchor1 <- anchors(diff_loop_gi, type="first")
anchor2 <- anchors(diff_loop_gi, type="second")

pbs_ld_anchor1_overlap <- findOverlaps(pbs_ld_GR, anchor1)
pbs_ld_anchor2_overlap <- findOverlaps(pbs_ld_GR, anchor2)
fnf_ld_anchor1_overlap <- findOverlaps(fnf_ld_GR, anchor1)
fnf_ld_anchor2_overlap <- findOverlaps(fnf_ld_GR, anchor2)

#hic_linked_pbs <- pbs_ld_GR[unique(queryHits(pbs_ld_anchor1_overlap), queryHits(pbs_ld_anchor2_overlap)),"peak"]
#hic_linked_fnf <- fnf_ld_GR[unique(queryHits(fnf_ld_anchor1_overlap), queryHits(fnf_ld_anchor2_overlap)),"peak"]
# c(hic_linked_pbs$peak, hic_linked_fnf$peak) |> unique() |> length()
# ============================================================
# annotate_longrange_ld_loop_genes
# ============================================================
annotate_longrange_ld_loop_genes <- function(ld_GR, ld_anchor1_overlap, ld_anchor2_overlap,
                                             loop_gi, promoters_gr, de_genes_df,
                                             condition_name, peak_region_map,
                                             padj_thresh = 0.05) {
  # Internal helper: keep_gene() now also returns symbol and genomic info
  keep_gene <- function(gid) {
    gid_clean <- gsub("\\..*", "", gid)
    hit <- de_genes_df[gene_id_noVer == gid_clean]
    if (nrow(hit) == 0) return(NULL)
    hit[1, .(
      is_de = !is.na(padj) & padj < padj_thresh,
      de_log2fc = log2FoldChange,
      de_padj = padj,
      de_direction = ifelse(log2FoldChange > 0, "UP_in_FNF", "DOWN_in_FNF"),
      gene_symbol = SYMBOL,
      seqnames = seqnames,
      start = start,
      end = end,
      strand = strand
    )]
  }

  results <- list()

  add_hits <- function(hits_obj, ld_anchor_label, gene_anchor_label) {
    if (length(hits_obj) == 0) return(invisible(NULL))
    qh <- queryHits(hits_obj)
    sh <- subjectHits(hits_obj)

    for (i in seq_along(qh)) {
      ld_idx   <- qh[i]
      loop_idx <- sh[i]
      peak_name <- ld_GR$peak[ld_idx]

      # region type (promoter, gene_body, distal)
      region_type <- peak_region_map[peak == peak_name, region_type]
      if (length(region_type) == 0) region_type <- NA

      lead_snp <- ld_GR$lead_caqtl[ld_idx]
      ld_snp   <- ld_GR$ld_snp[ld_idx]
      r2       <- ld_GR$R2[ld_idx]

      # get opposite anchor (gene side)
      gene_anchor <- anchors(loop_gi, type = ifelse(grepl("1$", ld_anchor_label),
                                                    "second", "first"))[loop_idx]

      promoter_ov <- findOverlaps(gene_anchor, promoters_gr)
      if (length(promoter_ov) == 0) next
      gene_ids <- promoters_gr$gene_id[subjectHits(promoter_ov)]

      for (gid in gene_ids) {
        de_info <- keep_gene(gid)
        if (is.null(de_info)) next

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
          gene_id        = gid,
          gene_symbol    = de_info$gene_symbol,
          seqnames       = de_info$seqnames,
          start          = de_info$start,
          end            = de_info$end,
          strand         = de_info$strand,
          is_de          = de_info$is_de,
          de_log2fc      = de_info$de_log2fc,
          de_padj        = de_info$de_padj,
          de_direction   = de_info$de_direction
        )
      }
    }
  }

  add_hits(ld_anchor1_overlap, "anchor1", "anchor2")
  add_hits(ld_anchor2_overlap, "anchor2", "anchor1")

  if (length(results) == 0) return(data.table())
  out <- rbindlist(results)

  # keep only DE genes (padj < threshold)
  out <- out[is_de == TRUE & !is.na(de_padj)]

  # deduplicate
  setorder(out, -R2)
  out <- out[!duplicated(out[, .(caqtl_peak, gene_id, loop_index)]), ]

  return(out[])
}

# ============================================================
# Run for PBS and FNF
# ============================================================

cat("Annotating PBS...\n")
pbs_longrange <- annotate_longrange_ld_loop_genes(
  ld_GR=pbs_ld_GR,
  ld_anchor1_overlap=pbs_ld_anchor1_overlap,
  ld_anchor2_overlap=pbs_ld_anchor2_overlap,
  loop_gi=diff_loop_gi,
  promoters_gr=promoters,
  de_genes_df=de_genes_clean,
  condition_name="PBS",
  peak_region_map=pbs_peak_map,
  padj_thresh=0.05)

cat("Annotating FNF...\n")
fnf_longrange <- annotate_longrange_ld_loop_genes(
  ld_GR=fnf_ld_GR,
  ld_anchor1_overlap=fnf_ld_anchor1_overlap,
  ld_anchor2_overlap=fnf_ld_anchor2_overlap,
  loop_gi=diff_loop_gi,
  promoters_gr=promoters,
  de_genes_df=de_genes_clean,
  condition_name="FNF",
  peak_region_map=fnf_peak_map,
  padj_thresh=0.05)

longrange_all <- rbindlist(list(pbs_longrange, fnf_longrange), fill=TRUE)
longrange_all[, .N, by = .(condition, region_type)]
# check whether it's also response caQTL

pbs_longrange_merged <- response_caQTL_pbs %>%
  inner_join(pbs_longrange, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

# For FNF
fnf_longrange_merged <- response_caQTL_fnf %>%
  inner_join(fnf_longrange, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))

#shared

shared_longrange_merged <- shared_fnf %>%
  inner_join(longrange_all, by = c("snp" = "lead_snp", "peak" = "caqtl_peak"))


# ============================================================
# save results
# ============================================================


save(pbs_longrange, fnf_longrange, longrange_all,
     file=file.path(plot_data_dir, "caQTL_Promoter_longrange.RData"))
fwrite(longrange_all, file.path(plot_data_dir, "caQTL_Promoter_longrange.tsv"), sep="\t")

save(pbs_longrange_merged, fnf_longrange_merged, shared_longrange_merged,
     file=file.path(plot_data_dir, "sig_response_caqtl_longrange_merged.RData"))
