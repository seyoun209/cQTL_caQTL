# ============================================================
# libs
# ============================================================
library(data.table)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(ggplot2)
library(ggtext)
library(colorspace)
library(plotgardener, lib.loc="/users/s/e/seyoun/R/dev")
library(RColorBrewer)
yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")

base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
response_dir <- file.path(base_dir, "caQTL/data/response_qtl")
plot_dir <- file.path(base_dir, "caQTL/plots")
plot_data_dir <- file.path(base_dir, "caQTL/data/plot_data")
# ============================================================
# inputs expected to already exist in your session:
#   pbs_rasqual_sig, fnf_rasqual_sig  # significant caQTLs (must have 'peak', 'snp')
#   pbs_eqtl, fnf_eqtl                # eQTL tables (must have 'variantID','gene_symbol')
#   pbs_sQTL, fnf_sQTL                # sQTL tables (must have 'var_id','phe_id','SYMBOL')
#   txdb                               # TxDb (hg38 or your choice)
#   ld_dir_eqtl                        # dir with eQTL LD files (no chr subdir required)
#   ld_dir_sqtl                        # dir with sQTL LD files (can be same as eQTL LD)
# ============================================================

# -------------------------------
# params
# -------------------------------
R2_THRESH <- 0.7
PROM_UP   <- 2500
PROM_DN   <- 500

# ============================================================
# helpers: genomic regions
# ============================================================
make_caQTL_GRanges <- function(peak_vec) {
  pv <- unique(na.omit(peak_vec))
  sp <- stringr::str_split_fixed(pv, "_", 3)
  GRanges(seqnames = sp[,1], ranges = IRanges(as.numeric(sp[,2]), as.numeric(sp[,3])))
}

classify_region <- function(peak_vec, promoters_gr, exons_gr, introns_gr) {
  # dedupe on peaks and map back
  pv <- unique(na.omit(peak_vec))
  gr <- make_caQTL_GRanges(pv)
  is_prom   <- countOverlaps(gr, promoters_gr) > 0
  is_exon   <- countOverlaps(gr, exons_gr) > 0
  is_intron <- countOverlaps(gr, introns_gr) > 0

  reg_map <- data.table(peak = pv,
                        region_type = ifelse(is_prom, "promoter",
                                             ifelse(is_exon | is_intron, "gene_body", "distal")))
  # return per input peak (including duplicates if any) for easy merge
  merge(data.table(peak = peak_vec), reg_map, by = "peak", all.x = TRUE)
}

# ============================================================
# helpers: LD file lookup + ID normalization
# ============================================================
# find per-SNP LD file either as "<ld_dir>/<lead>.ld" or "<ld_dir>/<chr>/<lead>.ld"
ld_path <- function(snp, ld_dir) {
  p1 <- file.path(ld_dir, paste0(snp, ".ld"))
  if (file.exists(p1)) return(p1)
  chr <- strsplit(snp, ":")[[1]][1]
  p2 <- file.path(ld_dir, chr, paste0(snp, ".ld"))
  if (file.exists(p2)) return(p2)
  return(NA_character_)
}

# strip to "chr:pos"
strip_to_pos <- function(x) sub(":[ACGT]+:[ACGT]+$", "", x)

# ============================================================
# annotate eQTL via caQTL LD buddies (R² ≥ threshold)
# match_by = "exact" (chr:pos:ref:alt) or "position" (chr:pos only)
# ============================================================
annotate_eqtl_ld <- function(caPeaks, eqtl, ld_dir, r2_threshold = 0.7,
                             match_by = c("exact","position")) {
  match_by <- match.arg(match_by)
  caPeaks <- as.data.table(caPeaks)
  eqtl    <- as.data.table(eqtl)

  if (!all(c("snp","peak") %in% names(caPeaks))) stop("caPeaks must have 'snp' and 'peak'")
  if (!all(c("variantID","gene_symbol") %in% names(eqtl))) stop("eqtl must have 'variantID','gene_symbol'")

  # preprocessed keys
  if (match_by == "position") {
    eqtl[, variant_key := strip_to_pos(variantID)]
  } else {
    eqtl[, variant_key := variantID]
  }

  caPeaks[, `:=`(eqtl_annotation = NA_character_, eqtl_match_count = 0L)]

  pb <- txtProgressBar(min = 0, max = nrow(caPeaks), style = 3)
  for (i in seq_len(nrow(caPeaks))) {
    snp <- caPeaks$snp[i]
    fp  <- ld_path(snp, ld_dir)
    if (is.na(fp)) { if (i %% 200 == 0) setTxtProgressBar(pb, i); next }

    ld <- fread(fp, showProgress = FALSE)
    if (!"R2" %in% names(ld) || !"SNP_B" %in% names(ld)) { if (i %% 200 == 0) setTxtProgressBar(pb, i); next }

    if (match_by == "position") {
      ld[, key := strip_to_pos(SNP_B)]
    } else {
      ld[, key := SNP_B]
    }

    ld_hi <- unique(ld[R2 >= r2_threshold, .(key, SNP_B, R2)])
    if (nrow(ld_hi) == 0) { if (i %% 200 == 0) setTxtProgressBar(pb, i); next }

    hits <- merge(eqtl[, .(variant_key, gene_symbol)], ld_hi, by.x = "variant_key", by.y = "key")
    if (nrow(hits) > 0) {
      caPeaks$eqtl_match_count[i] <- nrow(hits)
      caPeaks$eqtl_annotation[i]  <- paste0(unique(paste0(hits$SNP_B, ":", hits$gene_symbol,
                                                          " (R2=", round(hits$R2, 3), ")")),
                                            collapse = ";")
    }
    if (i %% 200 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  caPeaks[]
}

# ============================================================
# annotate sQTL via caQTL LD buddies (R² ≥ threshold)
# match_by = "exact" or "position" (same semantics)
# ============================================================
annotate_sqtl_ld <- function(caPeaks, sqtl, ld_dir, r2_threshold = 0.7,
                             match_by = c("exact","position")) {
  match_by <- match.arg(match_by)
  caPeaks <- as.data.table(caPeaks)
  sqtl    <- as.data.table(sqtl)

  if (!all(c("snp","peak") %in% names(caPeaks))) stop("caPeaks must have 'snp' and 'peak'")
  if (!all(c("var_id","phe_id","SYMBOL") %in% names(sqtl))) stop("sqtl must have 'var_id','phe_id','SYMBOL'")

  if (match_by == "position") {
    sqtl[, var_key := strip_to_pos(var_id)]
  } else {
    sqtl[, var_key := var_id]
  }

  caPeaks[, `:=`(sqtl_annotation = NA_character_, sqtl_match_count = 0L)]

  pb <- txtProgressBar(min = 0, max = nrow(caPeaks), style = 3)
  for (i in seq_len(nrow(caPeaks))) {
    snp <- caPeaks$snp[i]
    fp  <- ld_path(snp, ld_dir)
    if (is.na(fp)) { if (i %% 200 == 0) setTxtProgressBar(pb, i); next }

    ld <- fread(fp, showProgress = FALSE)
    if (!"R2" %in% names(ld) || !"SNP_B" %in% names(ld)) { if (i %% 200 == 0) setTxtProgressBar(pb, i); next }

    if (match_by == "position") {
      ld[, key := strip_to_pos(SNP_B)]
    } else {
      ld[, key := SNP_B]
    }

    ld_hi <- unique(ld[R2 >= r2_threshold, .(key, SNP_B, R2)])
    if (nrow(ld_hi) == 0) { if (i %% 200 == 0) setTxtProgressBar(pb, i); next }

    hits <- merge(sqtl[, .(var_key, phe_id, SYMBOL)], ld_hi, by.x = "var_key", by.y = "key")
    if (nrow(hits) > 0) {
      caPeaks$sqtl_match_count[i] <- nrow(hits)
      caPeaks$sqtl_annotation[i]  <- paste0(unique(paste0(hits$SNP_B, ":", hits$phe_id, ":", hits$SYMBOL,
                                                          " (R2=", round(hits$R2, 3), ")")),
                                            collapse = ";")
    }
    if (i %% 200 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  caPeaks[]
}

# ============================================================
# region-wise summaries (separate for PBS/FNF)
# ============================================================
regionwise_summary <- function(dt, type = c("eqtl","sqtl")) {
  type <- match.arg(type)
  count_col <- paste0(type, "_match_count")
  dt[, .(
    N   = .N,
    Yes = sum(get(count_col) > 0),
    No  = sum(get(count_col) == 0)
  ), by = region_type][order(match(region_type, c("promoter","gene_body","distal")))]
}

prep_plot_df <- function(sum_dt, dataset_label) {
  sum_dt[, dataset := dataset_label]
  long <- melt(sum_dt, id.vars = c("region_type","dataset","N"),
               variable.name = "annot", value.name = "count")
  long[, percentage := count / pmax(1, sum(count)), by = .(region_type, dataset)]
  long[, legend := paste(dataset, annot, sep = "_")]
  long[, label  := paste0(sprintf("%.1f%%", 100*percentage), "\nn=", count)]
  long
}

# plot_region_bars <- function(df, title_lab, palette_base) {
#   dataset <- unique(df$dataset)
#   cols <- setNames(
#     c(lighten(palette_base, 0.3), palette_base),
#     c(paste0(dataset, "_No"), paste0(dataset, "_Yes"))
#   ) 
#   ggplot(df, aes(x = region_type, y = 100*percentage, fill = legend)) +
#     geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.75) +
#     geom_text(aes(label = label),
#               position = position_stack(vjust = 0.5, reverse = TRUE),
#               size = 3, color = "black", lineheight = 0.9) +
#     scale_fill_manual(values = cols) +
#     labs(x = NULL, y = "% of caQTLs", title = title_lab) +
#     theme_minimal(base_family = "Helvetica") +
#     theme(
#       axis.text.x  = element_text(size = 9, color = "black"),
#       axis.text.y  = element_text(size = 9, color = "black"),
#       axis.title.y = element_text(size = 9, margin = margin(r = 5)),
#      panel.grid   = element_blank(),
#       legend.position = "none",
#       plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
#     )
# }


plot_region_bars <- function(df, title_lab, palette_base) {
  color_map <- c(
    "PBS_No"  = lighten("#2057A7", 0.3),
    "PBS_Yes" = "#2057A7",
    "FNF_No"  = lighten("#F2BC40", 0.3),
    "FNF_Yes" = "#F2BC40"
  )
  ggplot(df, aes(x = region_type, y = count, fill = legend, group = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = count),
              position = position_dodge(width = 0.8),
              size = 2, color = "black", lineheight = 0.9) +
    scale_fill_manual(values = color_map) +
    labs(x = NULL, y = "Number of caQTLs", title = title_lab) +
    theme_minimal(base_family = "Helvetica") +
    theme(plot.title = element_text(hjust = 0.5, size = 6),
          axis.text.x = element_text(size = 4),
          axis.text.y = element_text(size = 4),
          axis.title.y = element_text(size = 6),
          axis.line.y = element_line(size = 0.2, color = "black"),
          title = element_text(size=6),
          strip.placement = "outside",
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid   = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.x.bottom = element_text(size = 6),
          legend.key.size = unit(0.1, "cm"),
          legend.title = element_text(size = 4),
          legend.text  = element_text(size = 4))
}

# ============================================================
# 0) genomic regions
# ============================================================
txdb_genes   <- genes(txdb)
promoters_gr <- promoters(txdb_genes, upstream = PROM_UP, downstream = PROM_DN)
exons_gr     <- exons(txdb)
introns_gr   <- unlist(intronsByTranscript(txdb))

pbs_peak_map <- classify_region(pbs_rasqual_sig$peak, promoters_gr, exons_gr, introns_gr)
fnf_peak_map <- classify_region(fnf_rasqual_sig$peak, promoters_gr, exons_gr, introns_gr)

# ============================================================
# 1) eQTL: annotate & summarize (PBS / FNF separately)
#     choose match_by = "exact" or "position" (if alleles differ between resources)
# ============================================================
ld_dir_eqtl <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05"

pbs_eqtl_ld <- annotate_eqtl_ld(pbs_rasqual_sig, pbs_eqtl, ld_dir_eqtl,
                                r2_threshold = R2_THRESH, match_by = "exact")
fnf_eqtl_ld <- annotate_eqtl_ld(fnf_rasqual_sig, fnf_eqtl, ld_dir_eqtl,
                                r2_threshold = R2_THRESH, match_by = "exact")


save(pbs_peak_map, fnf_peak_map, pbs_eqtl_ld, fnf_eqtl_ld, file = file.path(plot_data_dir, "caQTL_eQTL_LD_annotated.Rdata"))
load(file.path(plot_data_dir, "caQTL_eQTL_LD_annotated.Rdata"))

pbs_eqtl_ld <- merge(pbs_eqtl_ld, pbs_peak_map, by = "peak", all.x = TRUE)
fnf_eqtl_ld <- merge(fnf_eqtl_ld, fnf_peak_map, by = "peak", all.x = TRUE)

pbs_eqtl_sum <- regionwise_summary(pbs_eqtl_ld, "eqtl")
fnf_eqtl_sum <- regionwise_summary(fnf_eqtl_ld, "eqtl")

caqtl_eqtl_pbs_pairID <- pbs_eqtl_ld |> filter(eqtl_match_count != 0 ) |> dplyr::select("Peak_SNP")
caqtl_eqtl_fnf_pairID <- fnf_eqtl_ld |> filter(eqtl_match_count != 0 ) |> dplyr::select("Peak_SNP")
rbind(caqtl_eqtl_pbs_pairID, caqtl_eqtl_fnf_pairID) |> unique() |> nrow()

pbs_eqtl_df <- prep_plot_df(pbs_eqtl_sum, "PBS")
fnf_eqtl_df <- prep_plot_df(fnf_eqtl_sum, "FNF")

# p_eqtl_pbs <- plot_region_bars(pbs_eqtl_df, "eQTL overlap (PBS)", "#2057A7")
# p_eqtl_fnf <- plot_region_bars(fnf_eqtl_df, "eQTL overlap (FNF)", "#F2BC40")

pbs_eqtl_df_yes <- pbs_eqtl_df[annot == "Yes"]
fnf_eqtl_df_yes <- fnf_eqtl_df[annot == "Yes"]
combined_eqtl_df_yes <- rbind(pbs_eqtl_df_yes, fnf_eqtl_df_yes)

combined_eqtl_df <- rbind(pbs_eqtl_df, fnf_eqtl_df)
p_eqtl_overlaps <- plot_region_bars(combined_eqtl_df_yes, "eQTL overlap (PBS & FNF)", "#2057A7")

save(p_eqtl_overlaps, file = file.path(plot_data_dir, "eQTL_overlaps.rds"))

# ============================================================
# 2) sQTL: annotate & summarize (PBS / FNF separately)
# ============================================================
ld_dir_sqtl <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05"

pbs_sqtl_ld <- annotate_sqtl_ld(pbs_rasqual_sig, pbs_sQTL, ld_dir_sqtl,
                                r2_threshold = R2_THRESH, match_by = "exact")
fnf_sqtl_ld <- annotate_sqtl_ld(fnf_rasqual_sig, fnf_sQTL, ld_dir_sqtl,
                                r2_threshold = R2_THRESH, match_by = "exact")

save(pbs_sqtl_ld, fnf_sqtl_ld, file = file.path(plot_data_dir, "caQTL_sQTL_LD_annotated.Rdata"))
load(file.path(plot_data_dir, "caQTL_sQTL_LD_annotated.Rdata"))

pbs_sqtl_ld <- merge(pbs_sqtl_ld, pbs_peak_map, by = "peak", all.x = TRUE)
fnf_sqtl_ld <- merge(fnf_sqtl_ld, fnf_peak_map, by = "peak", all.x = TRUE)

caqtl_sqtl_pbs_pairID <- pbs_sqtl_ld |> filter(sqtl_match_count != 0 ) |> dplyr::select("Peak_SNP")
caqtl_sqtl_fnf_pairID <- fnf_sqtl_ld |> filter(sqtl_match_count != 0 ) |> dplyr::select("Peak_SNP")
rbind(caqtl_sqtl_pbs_pairID, caqtl_sqtl_fnf_pairID) |> unique() |> nrow()

pbs_sqtl_sum <- regionwise_summary(pbs_sqtl_ld, "sqtl")
fnf_sqtl_sum <- regionwise_summary(fnf_sqtl_ld, "sqtl")

pbs_sqtl_df <- prep_plot_df(pbs_sqtl_sum, "PBS")
fnf_sqtl_df <- prep_plot_df(fnf_sqtl_sum, "FNF")

# p_sqtl_pbs <- plot_region_bars(pbs_sqtl_df, "sQTL overlap (PBS)", "#2057A7")
# p_sqtl_fnf <- plot_region_bars(fnf_sqtl_df, "sQTL overlap (FNF)", "#F2BC40")


pbs_sqtl_df_yes <- pbs_sqtl_df[annot == "Yes"]
fnf_sqtl_df_yes <- fnf_sqtl_df[annot == "Yes"]
combined_sqtl_df_yes <- rbind(pbs_sqtl_df_yes, fnf_sqtl_df_yes)

combined_sqtl_df <- rbind(pbs_sqtl_df, fnf_sqtl_df)

save(combined_eqtl_df, combined_sqtl_df, combined_eqtl_df_yes, combined_sqtl_df_yes, file = file.path(plot_data_dir, "combined_eqtl_sqtl_df.Rdata"))
load(file.path(plot_data_dir, "combined_eqtl_sqtl_df.Rdata"))
p_sqtl_overlaps <- plot_region_bars(combined_sqtl_df_yes, "sQTL overlap (PBS & FNF)", "#2057A7")
save(p_sqtl_overlaps, file = file.path(plot_data_dir, "sQTL_overlaps.rds"))

# visualize 

# ============================================================
PBS_hic <- "/proj/phanstiel_lab/Data/processed/CQTL/hic_caQTL/dietJuicer/output/PBS/PBS_inter_30.hic"
FNF_hic <- "/proj/phanstiel_lab/Data/processed/CQTL/hic_caQTL/dietJuicer/output/FNF/FNF_inter_30.hic"

pbs_atac <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/wasp/signals/merged/PBS_merged.bw"
fnf_atac <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/wasp/signals/merged/FNF_merged.bw"

pbs_rna <- "/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/signals/merged_norm/CTL_norm.bw"
fnf_rna <- "/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/signals/merged_norm/FNF_norm.bw"


# ============================================================
# Function

plot_caqtl_locus <- function(chrom, start, end, gene_label = NULL,
                             out_prefix = NULL, highlight_color = "grey40") {

  pageCreate(width = 6, height = 7, showGuides = FALSE)

  #  Hi-C parameters
  p <- pgParams(
    assembly = "hg38",
    chrom = chrom,
    chromstart = start - 1e5,
    chromend = end + 1e5,
    resolution = 10e3,
    norm = "SCALE",
    x = 0.6, y = 0.6, width = 4.8, height = 2
  )

  # ======================
  # PBS Hi-C
  dat <- readHic(file = PBS_hic, params = p)
  Z <- quantile(dat$counts[dat$counts > 0], .9)
  PBS_rect <- plotHicRectangle(params = p, data = PBS_hic, zrange = c(0, Z))
  plotText("PBS Hi-C", x = PBS_rect$x, y = PBS_rect$y, just = c("left","top"),
           fontcolor = "#2057A7", fontsize = 8, fontface = "bold")

  # FNF Hi-C
  dat <- readHic(file = FNF_hic, params = p)
  Z <- quantile(dat$counts[dat$counts > 0], .9)
  FNF_rect <- plotHicRectangle(params = p, data = FNF_hic, zrange = c(0, Z), y = "0.1b")
  plotText("FNF Hi-C", x = FNF_rect$x, y = FNF_rect$y, just = c("left","top"),
           fontcolor = "#F2BC40", fontsize = 8, fontface = "bold")

  # ======================
  # ATAC
  plotSignal(params = p, data = PBS_atac, y = "0.1b", height = 0.25,
             linecolor = "#2057A7", fill = "#2057A7", scale = FALSE,
             label = "ATAC PBS")
  plotSignal(params = p, data = FNF_atac, y = "0.1b", height = 0.25,
             linecolor = "#F2BC40", fill = "#F2BC40", scale = FALSE,
             label = "ATAC FNF")

  # RNA
  plotSignal(params = p, data = PBS_rna, y = "0.1b", height = 0.25,
             linecolor = lighten("#2057A7", 0.3), fill = lighten("#2057A7", 0.3),
             scale = FALSE, label = "RNA PBS")
  plotSignal(params = p, data = FNF_rna, y = "0.1b", height = 0.25,
             linecolor = lighten("#F2BC40", 0.3), fill = lighten("#F2BC40", 0.3),
             scale = FALSE, label = "RNA FNF")

  # ======================
  # Gene model
  plotGenes(params = p, height = 0.5, y = "0.1b")
  plotGenomeLabel(params = p, length = p$width, y = "0b")

  # Highlight
  annoHighlight(PBS_rect, chrom = chrom,
                chromstart = start, chromend = end,
                y = PBS_rect$y, height = 6.4,
                fill = NA, col = highlight_color, lwd = 1)

  if (!is.null(gene_label)) {
    plotText(gene_label, x = 0.4, y = 6.3,
             just = c("left","bottom"), fontsize = 9, fontface = "italic")
  }

  if (!is.null(out_prefix)) {
    pageSave(file = file.path(plot_dir, paste0(out_prefix, "_caQTL_signal.pdf")))
  }
}




#example
# load data
load(file.path(plot_data_dir, "caQTL_Promoter_longrange.RData"))#pbs_longrange, fnf_longrange, longrange_all,
load(file.path(plot_data_dir, "sig_response_caqtl_longrange_merged.RData")) #pbs_longrange_merged, fnf_longrange_merged, shared_longrange_merged, 

pdf(file.path(plot_dir, "poster_fig_hic.pdf"), width = 6.75, height = 8)
target_fnf <- fnf_longrange_merged[order(-R2)][1, ]
#target_fnf[, .(peak, snp, R2, gene_id, de_direction)]

pageCreate(width = 6.5, height = 8, default.units = "inches",showGuides = FALSE)
plot_caqtl_multitrack_with_hic(
  peak_id = target_fnf$peak,
  snp_id = target_fnf$snp,
  gene_id = target_fnf$gene_symbol,
  gene_start = target_fnf$start,
  gene_end = target_fnf$end,
  primary_dataset = "FNF",
  add_manhattan = TRUE,
  x_start = 0.5,
  y_start = 0.5,
  width = 2.5,
  height = 1.4,
  zoom_range = 400000,
  add_atac_signal = TRUE,
  add_rna_signal = TRUE)

target_pbs <- pbs_longrange_merged[order(-R2)][1, ]

plot_caqtl_multitrack_with_hic(
  peak_id = target_pbs$peak,
  snp_id = target_pbs$snp,
  gene_id = target_pbs$gene_symbol,
  gene_start = target_pbs$start,
  gene_end = target_pbs$end,
  primary_dataset = "PBS",
  add_manhattan = TRUE,
  x_start = 3.75,
  y_start = 0.5,
  width = 2.5,
  height = 1.4,
  zoom_range = 350000,
  add_atac_signal = TRUE,
  add_rna_signal = TRUE)

#eqtl overlap barplot

load(file = file.path(plot_data_dir, "eQTL_overlaps.rds"))

plotGG(p_eqtl_overlaps, x = 0.5, y = 5.5, width = 3, height = 2)
load(file = file.path(plot_data_dir, "sQTL_overlaps.rds"))

plotGG(p_sqtl_overlaps, x = 3.75, y = 5.5, width = 3, height = 2)


dev.off()
