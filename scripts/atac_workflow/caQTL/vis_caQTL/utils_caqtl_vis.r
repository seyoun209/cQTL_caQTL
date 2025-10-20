# Load libraries ----------------------------------------------------------

# ATAC signal track-------------------------------------------------------------------
atac_signal_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/wasp/signals/merged"
pbs_atac <- file.path(atac_signal_dir, "PBS_merged.bw")
fnf_atac <- file.path(atac_signal_dir, "FNF_merged.bw")

# RNA signal track (n=100)------------------------------------------------------------

rna_signal_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/signals/merged_norm"
pbs_rna <- file.path(rna_signal_dir, "CTL_norm.bw")  # CTL = PBS
fnf_rna <- file.path(rna_signal_dir, "FNF_norm.bw")
# Functions----------------------------------------------------------------

#--------------------------------------------------------------------------
make_caQTL_GRanges <- function(peak_vec) {
  # Remove NA values and get unique peaks
  peaks <- unique(na.omit(peak_vec))

  # Split the peak string to extract chromosome, start, and end (assumes "chr_start_end")
  peak_split <- stringr::str_split_fixed(peaks, "_", 3)
  chr <- peak_split[, 1]
  start <- as.numeric(peak_split[, 2])
  end <- as.numeric(peak_split[, 3])

  # Create GRanges object
  gr <- GenomicRanges::GRanges(seqnames = chr,
                               ranges = IRanges::IRanges(start = start, end = end))
  return(gr)
}

#------------------------------------------------------------------------------
# extend version gr -----------------------------------------------------------
make_caQTL_GRanges_expanded <- function(peak_vec, window_kb = 25) {
  # Remove NA values and get unique peaks
  peaks <- unique(na.omit(peak_vec))
  
  # Split the peak string to extract chromosome, start, and end
  peak_split <- stringr::str_split_fixed(peaks, "_", 3)
  chr <- peak_split[, 1]
  start <- as.numeric(peak_split[, 2])
  end <- as.numeric(peak_split[, 3])
  
  # Calculate peak center
  peak_center <- floor((start + end) / 2)
  
  # Expand by ±window_kb
  window_bp <- window_kb * 1000
  expanded_start <- peak_center - window_bp
  expanded_end <- peak_center + window_bp
  
  # Make sure start is not less than 1
  expanded_start <- pmax(expanded_start, 1)
  
  # Create GRanges object with expanded regions
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = expanded_start, end = expanded_end)
  )
  return(gr)
}


#------------------------------------------------------------------------------
# getting the LD lookup caqtl (r2 > 0.7)---------------------------------------
get_caqtl_ld_variants <- function(caPeaks, ld_dir, r2_threshold = 0.5) {
  # caPeaks: data.table with caQTL info (must include column 'snp')
  # Returns: data.table with caQTL SNP, LD variants, and R2 values
  
  result_list <- list()
  
  for(i in seq_len(nrow(caPeaks))) {
    snp_id <- caPeaks$snp[i]
    ld_file <- file.path(ld_dir, paste0(snp_id, ".ld"))
    
    if(file.exists(ld_file)) {
      tryCatch({
        ld_data <- fread(ld_file)
        ld_hi <- ld_data[R2 > r2_threshold]
        
        if(nrow(ld_hi) > 0) {
          ld_hi$lead_snp <- snp_id
          ld_hi$peak <- caPeaks$peak[i]
          result_list[[i]] <- ld_hi
        }
      }, error = function(e) {
        message(paste("Error processing", snp_id, ":", e$message))
      })
    }
  }
  
  # Combine all results
  if(length(result_list) > 0) {
    return(rbindlist(result_list, fill = TRUE))
  } else {
    return(data.table())
  }
}


#---------------------------------------------------------------------------------

collapse_ld_variants <- function(ld_dt) {
  if (is.null(ld_dt) || nrow(ld_dt) == 0) return(ld_dt)
  ld_dt %>%
    group_by(peak, lead_snp) %>%
    arrange(desc(R2)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    as.data.table()
}

convert_ld_to_granges <- function(ld_dt) {
  if (is.null(ld_dt) || nrow(ld_dt) == 0) return(GRanges())
  snp_split <- stringr::str_split_fixed(ld_dt$SNP_B, ":", 4)
  GRanges(
    seqnames = snp_split[, 1],
    ranges   = IRanges(start = as.numeric(snp_split[, 2]), width = 1),
    ld_snp       = ld_dt$SNP_B,
    lead_caqtl   = ld_dt$lead_snp,
    peak         = ld_dt$peak,
    R2           = ld_dt$R2
  )
}

# This is the annotate whether it is promoter vs distal using the peak span
label_peak_region_type <- function(peaks_char, promoters_gr) {
  gr0 <- make_caQTL_GRanges(peaks_char)
  ov  <- countOverlaps(gr0, promoters_gr) > 0
  data.table(peak = unique(na.omit(peaks_char)),
             region_type = ifelse(ov, "promoter", "distal"))
}

#--------------------------------------------------------------------------------

annotate_longrange_ld_loop_genes <- function(ld_GR,
                                             ld_anchor1_overlap,
                                             ld_anchor2_overlap,
                                             loop_gi,
                                             promoters_gr,
                                             de_genes_df,
                                             condition_name,
                                             peak_region_map,
                                             padj_thresh = 0.05) {

  keep_gene <- function(gid) {
    gid_clean <- gsub("\\..*", "", gid)
    hit <- de_genes_df[gene_id_noVer == gid_clean]
    if (nrow(hit) == 0) return(NULL)
    hit[1, .(is_de = !is.na(padj) & padj < padj_thresh,
             de_log2fc = log2FoldChange, de_padj = padj,
             de_direction = ifelse(log2FoldChange > 0, "UP_in_FNF", "DOWN_in_FNF"))]
  }

  results <- list()

  # helper to add records (LD at anchor X, gene at the other)
  add_hits <- function(hits_obj, ld_anchor_label, gene_anchor_label) {
    if (length(hits_obj) == 0) return(invisible(NULL))

    qh <- queryHits(hits_obj)
    sh <- subjectHits(hits_obj)

    for (i in seq_along(qh)) {
      ld_idx   <- qh[i]
      loop_idx <- sh[i]

      peak_name <- ld_GR$peak[ld_idx]
      # enforce LONG-RANGE: caQTL peak must be distal (not promoter-overlapping)
      if (!peak_name %in% peak_region_map$peak) next
      if (peak_region_map[peak == peak_name, region_type] != "distal") next

      lead_snp <- ld_GR$lead_caqtl[ld_idx]
      ld_snp   <- ld_GR$ld_snp[ld_idx]
      r2       <- ld_GR$R2[ld_idx]

      # the gene-side anchor
      gene_anchor <- anchors(loop_gi, type = ifelse(grepl("1$", ld_anchor_label), "second", "first"))[loop_idx]
      promoter_ov <- findOverlaps(gene_anchor, promoters_gr)

      if (length(promoter_ov) == 0) next
      gene_ids <- promoters_gr$gene_id[subjectHits(promoter_ov)]

      for (gid in gene_ids) {
        de_info <- keep_gene(gid)
        if (is.null(de_info)) next

        results[[length(results) + 1]] <<- data.table(
          condition      = condition_name,
          caqtl_peak     = peak_name,
          lead_snp       = lead_snp,
          ld_variant     = ld_snp,
          R2             = r2,
          loop_index     = loop_idx,
          ld_at_anchor   = ld_anchor_label,
          gene_at_anchor = gene_anchor_label,
          gene_id        = gid,
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

  # keep only DE genes (padj < thresh)
  out <- out[is_de == TRUE & !is.na(de_padj)]

  # dedupe per (peak, gene, loop) by the strongest LD (highest R2)
  setorder(out, -R2)
  out <- out[!duplicated(out[, .(caqtl_peak, gene_id, loop_index)]), ]
  out[]
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


plot_caqtl_multitrack_with_hic <- function(
  peak_id, snp_id, gene_id,
  gene_start, gene_end,
  primary_dataset = "PBS",
  x_start = 0.5, y_start = 0.5,
  width = 6, height = 5,
  zoom_range = 200000,
  add_atac_signal = TRUE, add_rna_signal = TRUE,
  add_hic = TRUE, add_manhattan = TRUE
) {

  if (!primary_dataset %in% c("PBS", "FNF"))
    stop("primary_dataset must be PBS or FNF")

  highConf_resQtL <- if (primary_dataset == "PBS") response_caQTL_pbs else response_caQTL_fnf
  test_boxplotInfo <- highConf_resQtL |> filter(peak == peak_id, snp == snp_id)
  if (nrow(test_boxplotInfo) == 0)
    stop("Peak–SNP pair not found in high-confidence QTLs")

  # ---- region
  parts <- strsplit(peak_id, "_")[[1]]
  chrom <- parts[1]
  peak_start <- as.numeric(parts[2])  # CHANGED from 'start'
  peak_end <- as.numeric(parts[3])    # CHANGED from 'end'
  center <- round((peak_start + peak_end) / 2)
  minregion <- max(1, center - zoom_range)
  maxregion <- center + zoom_range

region_pg <- pgParams(
  assembly = "hg38",
  chrom = chrom,
  chromstart = minregion,
  chromend = maxregion,
  x = x_start,
  width = width
)

  # ---- load LD + QTL
ld_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld0"
ld_file <- file.path(ld_dir, paste0(snp_id, ".ld"))
ld_full <- if (file.exists(ld_file)) fread(ld_file) else data.frame()

if ("SNP_B" %in% colnames(ld_full)) {
  ld_calc <- ld_full |> dplyr::select(SNP_B, R2) |> dplyr::rename(snp = SNP_B)
} else if ("snp" %in% colnames(ld_full)) {
  ld_calc <- ld_full |> dplyr::select(snp, R2)
} else {
  ld_calc <- data.frame(snp = character(), R2 = numeric())
}

pbs_qtl_region <- pbs_rasqual |> filter(Chromosome == chrom, peak == peak_id)
fnf_qtl_region <- fnf_rasqual |> filter(Chromosome == chrom, peak == peak_id)
leftjoin_pbs <- inner_join(pbs_qtl_region, ld_calc, by = "snp")
leftjoin_fnf <- inner_join(fnf_qtl_region, ld_calc, by = "snp")

prep_locus <- function(data) {
  data |>
    mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) |>
    mutate(LDgrp = addNA(LDgrp)) |>
    separate(snp, into = c("chrom", "pos", "ref", "alt"), sep = ":", remove = FALSE) |>
    mutate(pos = as.numeric(pos)) |>
    dplyr::rename(p = "PValue") |>
    filter(!is.na(p)) |>
    dplyr::select(chrom, pos, p, snp, R2, LDgrp)|>
        mutate(
        LDgrp = factor(LDgrp, 
                        levels = c(NA, "(0.8,1]", "(0.6,0.8]", "(0.4,0.6]", 
                                "(0.2,0.4]", "(0,0.2]"), 
                        ordered = TRUE)
        ) |>
        filter(!is.na(LDgrp)) |>
        arrange(desc(LDgrp)) |>
        data.frame()
}


pbs_locus <- prep_locus(leftjoin_pbs)
fnf_locus <- prep_locus(leftjoin_fnf)
ylim_pg <- max(
  ceiling(max(-log10(pbs_locus$p), na.rm = TRUE)) + 2,
  ceiling(max(-log10(fnf_locus$p), na.rm = TRUE)) + 2,
  5
)

  # ---- layout
track_gap <- 0.125
current_y <- y_start
track_gap <- 0.125
track_height <- 0.4
manhattan_height <- 0.5

# 1️⃣ Manhattan ------------------------------------------------
if (add_manhattan) {
  for (cond in c("PBS", "FNF")) {
  locus_data <- if (cond == "PBS") pbs_locus else fnf_locus
  y_offset <- current_y
  manhattan <- plotManhattan(
    data = locus_data,
    params = region_pg,
    range = c(0, ylim_pg),
    fill = colorby("LDgrp", palette = colorRampPalette(c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"))),
    y = y_offset,
    height = manhattan_height,
    snpHighlights = data.frame(snp = snp_id, pch = 24, cex = 0.5, col = "black")
  )
  plotText(label = paste(cond, "caQTL"),
            x = x_start + 0.1, y = y_offset  - 0.05,
            fontsize = 7, fontfamily = "Helvetica", just = c("left", "top"))
  annoYaxis(plot = manhattan, at = seq(0, ylim_pg, 2), axisLine = TRUE, fontsize = 6)
  current_y <- current_y + manhattan_height + track_gap 
  }
}

# Add y-axis label
plotText(
label = "-log10(p-value)",
x = x_start - 0.31, 
y = y_start + height / 2 ,
rot = 90, 
fontsize = 8, 
just = "center",
default.units = "inches", 
fontfamily = "Helvetica"
)

plotText(
label = snp_id,
x = x_start + width / 2, 
y = y_start - 0.1,
just = c("center", "top"),
fontfamily = "Helvetica", 
fontsize = 8
)
  # 2️⃣ Hi-C -------------------------------------------------------
  if (add_hic) {
    hic_params <- pgParams(
      assembly = "hg38", chrom = chrom,
      chromstart = minregion, chromend = maxregion,
      resolution = 10e3, norm = "SCALE",
      x = x_start, width = width
    )
    for (cond in c("PBS", "FNF")) {
      hic_file <- if (cond == "PBS") PBS_hic else FNF_hic
      col <- if (cond == "PBS") "#2057A7" else "#F2BC40"
      dat <- readHic(file = hic_file, params = hic_params)
      Z <- quantile(dat$counts[dat$counts > 0], .9)
      plotHicRectangle(params = hic_params, data = hic_file,
                       zrange = c(0, Z), y = current_y, height = 0.8)
      plotText(paste0(cond, " Hi-C"), x = x_start + 0.1, y = current_y - 0.05,
               fontsize = 7, fontcolor = col, fontface = "bold")
      current_y <- current_y + 0.8 + track_gap
    }
  }

  # 3️⃣ ATAC -------------------------------------------------------
  if (add_atac_signal) {
  ATAC_signals <- plotMultiSignal(
      data = list(pbs_atac, fnf_atac),  # Now uses file paths defined above
      params = region_pg,
      y = current_y - 0.1,
      height = track_height,
      linecolor = c(yl_gn_bu[8], yl_gn_bu[8]),
      fill = c(yl_gn_bu[8], yl_gn_bu[8]),
      default.units = "inches",
      gapdistance = 0.04
    )
    plotText(label = "ATAC", x = x_start - 0.3, y = current_y + track_height / 2 - 0.02,
             rot = 90, fontsize = 8, just = "center", fontfamily = "Helvetica", fontcolor = yl_gn_bu[8])
    plotText(label = "PBS", x = x_start - 0.05, y = current_y - 0.02,
             fontsize = 6, fontfamily = "Helvetica", fontcolor = yl_gn_bu[8])
    plotText(label = "FN-f", x = x_start - 0.05, y = current_y + track_height / 2 - 0.02,
             fontsize = 6, fontfamily = "Helvetica", fontcolor = yl_gn_bu[8])
    current_y <- current_y + track_height
  }
  annoHighlight(
  plot= ATAC_signals[[1]],
  chrom = chrom, 
  chromstart = peak_start-1000,
  chromend = peak_end+1000,
  y = current_y-track_height-0.1 , height = 1.5,
  default.units = "inches", just = c("left", "top")
)


# -----------------------------------------------------------
# 4️⃣ RNA track with DE-gene highlights
# -----------------------------------------------------------
if (add_rna_signal) {
RNA_signals <- plotMultiSignal(
  data = list(pbs_rna, fnf_rna),
  params = region_pg,
  y = current_y,
  height = track_height,
  linecolor = c(yl_gn_bu[6], yl_gn_bu[6]),
  fill = c(yl_gn_bu[6], yl_gn_bu[6]),
  default.units = "inches",
  gapdistance = 0.04
)
plotText(label = "RNA", x = x_start - 0.3, y = current_y + track_height / 2 + 0.02,
          rot = 90, fontsize = 8, just = "center", fontfamily = "Helvetica", fontcolor = yl_gn_bu[6])
plotText(label= "PBS", x = x_start-0.05, y = current_y + 0.04,
      fontsize=6, fontfamily= "Helvetica", fontcolor = yl_gn_bu[6])
plotText(label= "FN-f", x = x_start-0.05, y = current_y + track_height/2 + 0.04,
      fontsize=6, fontfamily= "Helvetica", fontcolor = yl_gn_bu[6])
current_y <- current_y + track_height
}

annoHighlight(
plot= RNA_signals[[1]],
chrom = chrom,
chromstart = gene_start-1000,
chromend = gene_start+1000,
y = current_y - (track_height *2) - 0.1, height = 1.5,
default.units = "inches", just = c("left", "top")
)


# -----------------------------------------------------------
# 5️⃣ Gene annotation track (with gene highlights)
# -----------------------------------------------------------
gene_highlights <- data.frame("gene" = gene_id, color = "#225EA8")


gene_plot <- plotGenes(
  params = region_pg,
  y = current_y,
  height = 0.4,
  fontsize = 6,
  geneHighlights = if (nrow(gene_highlights) > 0) gene_highlights else NULL,
  geneBackground = "grey"
)

annoGenomeLabel(
  plot = gene_plot,
  params = region_pg,
  fontsize = 6,
  y = current_y + 0.4
)
}




#======================================================================
#eqtl -caqtl -hic
#======================================================================

plot_caqtl_multitrack_with_hic_eqtl <- function(
  peak_id, snp_id, 
  gene_id,         # ENSG ID for eQTL filtering
  gene_symbol,     # Gene symbol for display/highlighting
  gene_start, gene_end,
  primary_dataset = "PBS",
  eqtl_dataset = NULL,
  x_start = 0.5, y_start = 0.5,
  width = 6, height = 5,
  zoom_range = 200000,
  add_atac_signal = TRUE, 
  add_rna_signal = TRUE,
  add_hic = TRUE, 
  add_manhattan = TRUE,
  add_eqtl = FALSE
) {

  if (!primary_dataset %in% c("PBS", "FNF"))
    stop("primary_dataset must be PBS or FNF")
  
  # Default eqtl_dataset to primary_dataset if add_eqtl is TRUE
  if (add_eqtl && is.null(eqtl_dataset)) {
    eqtl_dataset <- primary_dataset
  }

  highConf_resQtL <- if (primary_dataset == "PBS") pbs_caqtl_hic_eqtl_response else fnf_caqtl_hic_eqtl_response
  test_boxplotInfo <- highConf_resQtL |> filter(peak == peak_id, snp == snp_id)
  if (nrow(test_boxplotInfo) == 0)
    stop("Peak–SNP pair not found in high-confidence QTLs")

  # ---- region
  parts <- strsplit(peak_id, "_")[[1]]
  chrom <- parts[1]
  peak_start <- as.numeric(parts[2])
  peak_end <- as.numeric(parts[3])
  center <- round((peak_start + peak_end) / 2)
  minregion <- max(1, center - zoom_range)
  maxregion <- center + zoom_range

  region_pg <- pgParams(
    assembly = "hg38",
    chrom = chrom,
    chromstart = minregion,
    chromend = maxregion,
    x = x_start,
    width = width
  )

  # ---- load LD + QTL (caQTL data)
  ld_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld0"
  ld_file <- file.path(ld_dir, paste0(snp_id, ".ld"))
  ld_full <- if (file.exists(ld_file)) fread(ld_file) else data.frame()

  if ("SNP_B" %in% colnames(ld_full)) {
    ld_calc <- ld_full |> dplyr::select(SNP_B, R2) |> dplyr::rename(snp = SNP_B)
  } else if ("snp" %in% colnames(ld_full)) {
    ld_calc <- ld_full |> dplyr::select(snp, R2)
  } else {
    ld_calc <- data.frame(snp = character(), R2 = numeric())
  }

  pbs_qtl_region <- pbs_rasqual |> filter(Chromosome == chrom, peak == peak_id)
  fnf_qtl_region <- fnf_rasqual |> filter(Chromosome == chrom, peak == peak_id)
  leftjoin_pbs <- inner_join(pbs_qtl_region, ld_calc, by = "snp")
  leftjoin_fnf <- inner_join(fnf_qtl_region, ld_calc, by = "snp")

  prep_locus <- function(data) {
    data |>
      mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) |>
      mutate(LDgrp = addNA(LDgrp)) |>
      separate(snp, into = c("chrom", "pos", "ref", "alt"), sep = ":", remove = FALSE) |>
      mutate(pos = as.numeric(pos)) |>
      dplyr::rename(p = "PValue") |>
      filter(!is.na(p)) |>
      dplyr::select(chrom, pos, p, snp, R2, LDgrp) |>
      mutate(
        LDgrp = factor(LDgrp, 
                      levels = c(NA, "(0.8,1]", "(0.6,0.8]", "(0.4,0.6]", 
                                "(0.2,0.4]", "(0,0.2]"), 
                      ordered = TRUE)
      ) |>
      filter(!is.na(LDgrp)) |>
      arrange(desc(LDgrp)) |>
      data.frame()
  }

  pbs_locus <- prep_locus(leftjoin_pbs)
  fnf_locus <- prep_locus(leftjoin_fnf)
  ylim_pg <- max(
    ceiling(max(-log10(pbs_locus$p), na.rm = TRUE)) + 2,
    ceiling(max(-log10(fnf_locus$p), na.rm = TRUE)) + 2,
    5
  )

  # ---- layout
  track_gap <- 0.125
  current_y <- y_start
  track_height <- 0.4
  manhattan_height <- 0.5

  # ============================================================
  # 1️⃣ caQTL Manhattan plots
  # ============================================================
  if (add_manhattan) {
    for (cond in c("PBS", "FNF")) {
      locus_data <- if (cond == "PBS") pbs_locus else fnf_locus
      y_offset <- current_y
      manhattan <- plotManhattan(
        data = locus_data,
        params = region_pg,
        range = c(0, ylim_pg),
        fill = colorby("LDgrp", palette = colorRampPalette(c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"))),
        y = y_offset,
        height = manhattan_height,
        snpHighlights = data.frame(snp = snp_id, pch = 24, cex = 0.5, col = "black")
      )
      plotText(label = paste(cond, "caQTL"),
              x = x_start + 0.1, y = y_offset - 0.05,
              fontsize = 7, fontfamily = "Helvetica", just = c("left", "top"))
      annoYaxis(plot = manhattan, at = seq(0, ylim_pg, 2), axisLine = TRUE, fontsize = 6)
      current_y <- current_y + manhattan_height + track_gap 
    }
    
    # Add y-axis label for caQTL
    plotText(
      label = "-log10(p)",
      x = x_start - 0.25, 
      y = y_start + manhattan_height,
      rot = 90, 
      fontsize = 7, 
      just = "center",
      default.units = "inches", 
      fontfamily = "Helvetica"
    )
    
    plotText(
      label = paste0("caQTL SNP: ", snp_id),
      x = x_start + width / 2, 
      y = y_start - 0.1,
      just = c("center", "top"),
      fontfamily = "Helvetica", 
      fontsize = 6
    )
  }

  # ============================================================
  # 1.5️⃣ eQTL Manhattan plots
  # ============================================================
  if (add_eqtl) {
  eqtl_nom_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/eqtl/qtl_nom"
  chr_tag <- sub("^chr", "", chrom)
  
  # Determine which conditions to plot
  eqtl_conditions <- if (eqtl_dataset == "both") c("PBS", "FNF") else eqtl_dataset
  
  # Clean gene_id for matching (remove version if present)
  gene_id_clean <- sub("\\..*", "", gene_id)
  
  for (cond in eqtl_conditions) {
    # Find eQTL file
    base <- if (cond == "PBS") "CTL_PEER_k20_genoPC" else "FNF_PEER_k22_genoPC"
    eqtl_file <- file.path(eqtl_nom_dir, sprintf("%s_allSignals_nom1Mb_MAFs_chr%s.csv", base, chr_tag))
    if (!file.exists(eqtl_file))
      eqtl_file <- file.path(eqtl_nom_dir, sprintf("%s_nom1Mb_chr%s.csv", base, chr_tag))
    if (!file.exists(eqtl_file)) {
      warning("No eQTL file found for ", cond, " ", chrom, ". Skipping eQTL track.")
      next
    }
    
    message("Reading eQTL data: ", basename(eqtl_file))
    eqtl_df <- fread(eqtl_file)
    
    # Filter to gene of interest using gene_id (ENSG)
    eqtl_gene <- eqtl_df %>%
      mutate(gene_id_clean = sub("\\..*", "", gene_id)) %>%
      filter(gene_id_clean == !!gene_id_clean) %>%
      mutate(
        chrom = variant_chr,
        pos = as.numeric(variant_start),
        p = as.numeric(nom_pval),
        R2 = as.numeric(r_squared),
      ) %>%
      filter(chrom == !!chrom & pos >= minregion & pos <= maxregion)
    
    if (nrow(eqtl_gene) == 0) {
      warning("No eQTL data for gene ", gene_id, " (", gene_symbol, ") in ", cond, " within region")
      next
    }
    
    message("Found ", nrow(eqtl_gene), " eQTL variants for ", gene_symbol, " (", gene_id, ")")
    
    # Prepare LD groups and MATCH caQTL locus format exactly
    eqtl_locus <- eqtl_gene %>%
      mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) %>%
      mutate(LDgrp = addNA(LDgrp)) %>%
      mutate(
        LDgrp = factor(LDgrp, 
                      levels = c(NA, "(0.8,1]", "(0.6,0.8]", "(0.4,0.6]", 
                                "(0.2,0.4]", "(0,0.2]"), 
                      ordered = TRUE)
      ) %>%
      filter(!is.na(p)) %>%
      dplyr::rename("snp" = variantID) %>%
      # Select ONLY the columns needed for Manhattan plot, in the same order
      dplyr::select(chrom, pos, p, snp, R2, LDgrp) %>%
      arrange(desc(LDgrp)) %>%
      as.data.frame()
    
    message("Prepared ", nrow(eqtl_locus), " eQTL variants for plotting")
    
    # Set y-axis range
    ylim_eqtl <- max(
      ceiling(max(-log10(eqtl_locus$p), na.rm = TRUE)) + 2,
      5
    )
    
    # Find top eQTL variant for highlighting
    top_eqtl_idx <- which.min(eqtl_locus$p)
    if (length(top_eqtl_idx) > 0) {
      top_eqtl_var <- eqtl_locus$snp[top_eqtl_idx[1]]
    } else {
      top_eqtl_var <- NULL
    }
    
    # Plot eQTL Manhattan
    y_offset <- current_y
    
    # Create snpHighlights only if we have a top variant
    if (!is.null(top_eqtl_var) && !is.na(top_eqtl_var)) {
      snp_highlights <- data.frame(snp = top_eqtl_var, pch = 24, cex = 0.5, col = "black")
    } else {
      snp_highlights <- NULL
    }
    
    eqtl_manhattan <- plotManhattan(
      data = eqtl_locus,
      params = region_pg,
      range = c(0, ylim_eqtl),
      fill = colorby("LDgrp", palette = colorRampPalette(c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"))),
      y = y_offset,
      height = manhattan_height,
      snpHighlights = snp_highlights
    )
    
    plotText(label = paste(cond, "eQTL"),
            x = x_start + 0.1, y = y_offset - 0.05,
            fontsize = 7, fontfamily = "Helvetica", just = c("left", "top"))
    annoYaxis(plot = eqtl_manhattan, at = seq(0, ylim_eqtl, 2), axisLine = TRUE, fontsize = 6)
    current_y <- current_y + manhattan_height + track_gap
  }
  
  # Add gene label for eQTL (use gene_symbol for display)
  plotText(
    label = paste0("eQTL gene: ", gene_symbol),
    x = x_start + width / 2, 
    y = current_y - manhattan_height - track_gap - 0.05,
    just = c("center", "top"),
    fontfamily = "Helvetica", 
    fontsize = 6,
    fontface = "italic"
  )
}
  # ============================================================
  # 2️⃣ Hi-C
  # ============================================================
  if (add_hic) {
    hic_params <- pgParams(
      assembly = "hg38", chrom = chrom,
      chromstart = minregion, chromend = maxregion,
      resolution = 10e3, norm = "SCALE",
      x = x_start, width = width
    )
    for (cond in c("PBS", "FNF")) {
      hic_file <- if (cond == "PBS") PBS_hic else FNF_hic
      col <- if (cond == "PBS") "#2057A7" else "#F2BC40"
      dat <- readHic(file = hic_file, params = hic_params)
      Z <- quantile(dat$counts[dat$counts > 0], .9)
      plotHicRectangle(params = hic_params, data = hic_file,
                       zrange = c(0, Z), y = current_y, height = 0.8)
      plotText(paste0(cond, " Hi-C"), x = x_start + 0.1, y = current_y - 0.05,
               fontsize = 7, fontcolor = col, fontface = "bold")
      current_y <- current_y + 0.8 + track_gap
    }
  }

  # ============================================================
  # 3️⃣ ATAC
  # ============================================================
if (add_atac_signal) {
  ATAC_signals <- plotMultiSignal(
    data = list(pbs_atac, fnf_atac),
    params = region_pg,
    y = current_y - 0.1,
    height = track_height,
    linecolor = c(yl_gn_bu[8], yl_gn_bu[8]),
    fill = c(yl_gn_bu[8], yl_gn_bu[8]),
    default.units = "inches",
    gapdistance = 0.04
  )
  plotText(label = "ATAC", x = x_start - 0.2, y = current_y + track_height / 2 - 0.02,
           rot = 90, fontsize = 7, just = "center", fontfamily = "Helvetica", fontcolor = yl_gn_bu[8])
  plotText(label = "PBS", x = x_start - 0.05, y = current_y - 0.02,
           fontsize = 5, fontfamily = "Helvetica", fontcolor = yl_gn_bu[8])
  plotText(label = "FN-f", x = x_start - 0.05, y = current_y + track_height / 2 - 0.02,
           fontsize = 5, fontfamily = "Helvetica", fontcolor = yl_gn_bu[8])
  current_y <- current_y + track_height
  
  # Add highlight only if peak coordinates are valid
  if (!is.null(peak_start) && !is.null(peak_end) && 
      length(peak_start) > 0 && length(peak_end) > 0 &&
      !is.na(peak_start) && !is.na(peak_end)) {
    annoHighlight(
      plot = ATAC_signals[[1]],
      chrom = chrom,
      chromstart = peak_start - 1000,
      chromend = peak_end + 1000,
      y = current_y - track_height - 0.1, 
      height = 1.5,
      default.units = "inches", 
      just = c("left", "top")
    )
  }
}

  # ============================================================
  # 4️⃣ RNA track
  # ============================================================
 if (add_rna_signal) {
  RNA_signals <- plotMultiSignal(
    data = list(pbs_rna, fnf_rna),
    params = region_pg,
    y = current_y,
    height = track_height,
    linecolor = c(yl_gn_bu[6], yl_gn_bu[6]),
    fill = c(yl_gn_bu[6], yl_gn_bu[6]),
    default.units = "inches",
    gapdistance = 0.04
  )
  plotText(label = "RNA", x = x_start - 0.2, y = current_y + track_height / 2 + 0.02,
           rot = 90, fontsize = 7, just = "center", fontfamily = "Helvetica", fontcolor = yl_gn_bu[6])
  plotText(label = "PBS", x = x_start - 0.05, y = current_y + 0.04,
           fontsize = 5, fontfamily = "Helvetica", fontcolor = yl_gn_bu[6])
  plotText(label = "FN-f", x = x_start - 0.05, y = current_y + track_height / 2 + 0.04,
           fontsize = 5, fontfamily = "Helvetica", fontcolor = yl_gn_bu[6])
  current_y <- current_y + track_height
  
  # Add highlight only if gene coordinates are valid
  if (!is.null(gene_start) && !is.null(gene_end) && 
      length(gene_start) > 0 && length(gene_end) > 0 &&
      !is.na(gene_start) && !is.na(gene_end)) {
    annoHighlight(
      plot = RNA_signals[[1]],
      chrom = chrom,
      chromstart = gene_start - 1000,
      chromend = gene_end + 1000,
      y = current_y - (track_height * 2) - 0.1, 
      height = 1.5,
      default.units = "inches", 
      just = c("left", "top")
    )
  }
}
  # ============================================================
  # 5️⃣ Gene annotation track (use gene_symbol for highlighting)
  # ============================================================
  gene_highlights <- data.frame("gene" = gene_symbol, color = yl_gn_bu[8])

  gene_plot <- plotGenes(
    params = region_pg,
    y = current_y,
    height = 0.4,
    fontsize = 5,
    geneHighlights = if (nrow(gene_highlights) > 0) gene_highlights else NULL,
    geneBackground = "grey"
  )

  annoGenomeLabel(
    plot = gene_plot,
    params = region_pg,
    fontsize = 5,
    y = current_y + 0.4
  )
}
