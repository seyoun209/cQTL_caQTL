library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggtext)
library(lme4)
library(plotgardener, lib.loc="/users/s/e/seyoun/R/dev")
library(RColorBrewer)
yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")
# Load data ----------------------------------------------------------
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
response_dir <- file.path(base_dir, "caQTL/data/response_qtl")
load(file.path(response_dir, "caQTL_prepdata.RData")) #vst_counts, peak_info,  all_rasqual_qtl_duplicated,all_caqtl_nodup,  all_caqtl, combined_geno_matrix, full_covar, meta_final
load(file.path(response_dir, "response_caQTL_PBS_highconf.RData"))
load(file.path(response_dir, "response_caQTL_FNF_highconf.RData"))

# nominal rds
pbs_rasqual <- readRDS(file.path(response_dir, "pbs_rasqual_nominal.rds"))
fnf_rasqual <- readRDS(file.path(response_dir, "fnf_rasqual_nominal.rds"))

# 100kb nominal rds
rasqual_dir <- file.path(base_dir, "caQTL/rasqual_output",paste0("combined_window_", 100,"kb"))
pbs_rasqual <- fread(file.path(rasqual_dir, "pbs_pc0_100kb_combined.txt"))
pbs_rasqual <-  pbs_rasqual |> mutate(Condition = "pbs") |>
  dplyr::rename(peak = Feature, snp = rs_ID) |>
  mutate(Peak_SNP = paste(peak, snp, sep = "_"))

fnf_rasqual <- fread(file.path(rasqual_dir, "fnf_pc0_100kb_combined.txt"))
fnf_rasqual <-   fnf_rasqual |> mutate(Condition = "fnf") |> 
    dplyr::rename(peak = Feature, snp = rs_ID) |>
    mutate(Peak_SNP = paste(peak, snp, sep = "_"))

pbs_rasqual_100k <- pbs_rasqual
fnf_rasqual_100k <- fnf_rasqual

#pbs_rasqual <- pbs_rasqual_100k
#fnf_rasqual <- fnf_rasqual_100k
meta_final <- meta_final |>
  mutate(
    Condition = ifelse(Condition == "CTL", 0,
                       ifelse(Condition == "FNF", 1, NA)),
    Donor = as.factor(Donor)
  ) |> 
  dplyr::select(-Sex)

# ATAC signal track-------------------------------------------------------------------
atac_signal_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/wasp/signals/merged"
pbs_atac <- file.path(atac_signal_dir, "PBS_merged.bw")
fnf_atac <- file.path(atac_signal_dir, "FNF_merged.bw")

# RNA signal track (n=100)------------------------------------------------------------

rna_signal_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/output/signals/merged_norm"
pbs_rna <- file.path(rna_signal_dir, "CTL_norm.bw")  # CTL = PBS
fnf_rna <- file.path(rna_signal_dir, "FNF_norm.bw")

# Functions ----------------------------------------------------------
# ============================================================================
# boxplot for caQTL
# ============================================================================
create_caqtl_boxplot <- function(peak_id, snp_id, highConf) {
    # Filter data for the specified peak-SNP pair
    test_boxplotInfo <- highConf |> 
        dplyr::filter(peak == peak_id, snp == snp_id)
    
    if (nrow(test_boxplotInfo) == 0) {
        stop("Peak-SNP pair not found in high-confidence response caQTLs")
    }
    
    # Extract normalized counts for this peak
    norm_counts_df <- vst_counts[peak_id, ] |>
        as.numeric() |>
        (\(x) data.frame(sampleID = colnames(vst_counts), vsd = x))()
    
    # Extract genotype data
    geno_row <- combined_geno_matrix |> filter(SNP == snp_id)
    if (nrow(geno_row) == 0) stop("SNP not found in genotype matrix")
    
    geno_samples <- setdiff(colnames(geno_row), 
                            c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"))
    matched_samples <- intersect(colnames(vst_counts), geno_samples)
    ordered_samples <- colnames(vst_counts)[colnames(vst_counts) %in% matched_samples]
    genotype_values <- as.numeric(unlist(geno_row[, ..ordered_samples, drop = FALSE]))
    geno_df <- data.frame(sampleID = ordered_samples, genotype = genotype_values)
    
    # Determine minor allele
    alleles <- unlist(strsplit(snp_id, ":"))
    ref_allele <- test_boxplotInfo$Ref
    alt_allele <- test_boxplotInfo$Alt
    
    # Get allele frequency to determine minor allele
    af <- test_boxplotInfo$Allele_Frequency
    minor_allele <- ifelse(af < 0.5, alt_allele, ref_allele)
    protective_allele <- ifelse(minor_allele == ref_allele, alt_allele, ref_allele)
    
    # Combine all data
    df <- Reduce(function(x, y) merge(x, y, by = "sampleID"),
                list(meta_final, norm_counts_df, geno_df, full_covar))
    df <- df[!is.na(df$genotype) & !is.na(df$vsd), ]
    
    # ---- Calculate RESIDUALS ----
    # Fit reduced model WITHOUT genotype (to see genotype effect in plot)
    pc_cols <- setdiff(colnames(full_covar), "sampleID")
    pc_rhs <- paste(pc_cols, collapse = " + ")
    
    reduced_formula <- as.formula(paste("vsd ~ Condition +", pc_rhs, "+ (1|Donor)"))
    reduced_model <- lmer(reduced_formula, data = df, REML = FALSE)
    
    # Extract residuals
    df$residual <- residuals(reduced_model)
    
    # Prepare plotting data
    meta_combined_all <- df %>%
        mutate(
        genotype_num = genotype,
        Condition = ifelse(Condition == 0, "PBS", "FN-f"),
        Condition = factor(Condition, levels = c("PBS", "FN-f")),
        genotype = case_when(
            genotype_num == 2 ~ paste(minor_allele, minor_allele, sep = "/"),
            genotype_num == 1 ~ paste(protective_allele, minor_allele, sep = "/"),
            genotype_num == 0 ~ paste(protective_allele, protective_allele, sep = "/")
        ),
        genotype = factor(genotype, levels = c(
            paste(protective_allele, protective_allele, sep = "/"),
            paste(protective_allele, minor_allele, sep = "/"),
            paste(minor_allele, minor_allele, sep = "/")
        ))
        )
  
  # Get effect sizes for text annotation

condition <- unique(test_boxplotInfo$Condition)

if (condition == "pbs") {
    pval_fnf <- fnf_rasqual |> filter(peak == peak_id, snp == snp_id) |> pull(PValue)
    pval_pbs <- test_boxplotInfo$RASQUAL_PValue_pbs[1]
    text_data <- data.frame(
        Condition = c("PBS", "FN-f"),
        beta = c(test_boxplotInfo$EffectSize_pbs[1], 
                 test_boxplotInfo$EffectSize_fnf[1]),
        pvalue = c(pval_pbs, pval_fnf)
    )
} else {
    pval_pbs <- pbs_rasqual |> filter(peak == peak_id, snp == snp_id) |> pull(PValue)
    pval_fnf <- test_boxplotInfo$RASQUAL_PValue_fnf[1]
    text_data <- data.frame(
        Condition = c("PBS", "FN-f"),
        beta = c(test_boxplotInfo$EffectSize_pbs[1], 
                 test_boxplotInfo$EffectSize_fnf[1]),
        pvalue = c(pval_pbs, pval_fnf)
    )
}
text_data <- data.frame(
    Condition = c("PBS", "FN-f"),
    beta = c(test_boxplotInfo$EffectSize_pbs[1], 
            test_boxplotInfo$EffectSize_fnf[1]),
    pvalue = c(pval_pbs, pval_fnf)
)

text_data <- text_data |> mutate(Condition = factor(Condition, levels = c("PBS", "FN-f")))
  
# Create plot
p <- ggplot(meta_combined_all, aes(x = genotype, y = residual, fill = Condition)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, alpha = 0.7) +
    scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
    geom_point(color = "grey40", 
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
               size = 0.25) +
    labs(x = "Genotype", y = "Chromatin accessibility (residuals)") +
    scale_fill_manual(values = c('#0067B9', '#FCCE52')) +
    scale_color_manual(values = c('#0067B9', '#FCCE52')) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
    strip.placement = "outside",
    axis.line.y = element_line(linewidth = 0.25),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 6, family = "Helvetica",margin = margin(r = -0.5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(0.5, "cm")) +
    facet_wrap(~ Condition, ncol = 2) +
    geom_text(data = text_data,aes(x = 2, y = Inf,
        label = ifelse(is.na(pvalue),
            sprintf("%s\nEffect: %.3f", Condition, beta),
            sprintf("%s\nEffect: %.3f\np: %.2e", Condition, beta, pvalue)),
        color = Condition),
        hjust = 0.5, vjust = 1, size = 2, family = "Helvetica", 
        lineheight = 0.75, check_overlap = TRUE
    )
  
return(p)
}

# ============================================================================
# manhattan plot for caQTL
# ============================================================================
plot_caqtl_multitrack <- function(peak_id, snp_id, primary_dataset = "PBS", x_start = 0.5, y_start = 0.5, width = 6, height = 5, zoom_range = 200000, add_atac_signal = TRUE, add_rna_signal = TRUE) {
# Validate input ----------------------------------------------------
if (!primary_dataset %in% c("PBS", "FNF")) {
  stop("Primary dataset must be either 'PBS' or 'FNF'")
}
highConf_resQtL <- if (primary_dataset == "PBS") pbs_highconf else fnf_highconf
test_boxplotInfo <- highConf_resQtL |> 
dplyr::filter(peak == peak_id, snp == snp_id)

if (nrow(test_boxplotInfo) == 0) stop("Peak–SNP pair not found in high-confidence set")

# Extract region info -----------------------------------------------
peak_parts <- unlist(strsplit(peak_id, "_"))
chrom <- peak_parts[1]
peak_start <- as.numeric(peak_parts[2])
peak_end <- as.numeric(peak_parts[3])
peak_center <- round((peak_start + peak_end) / 2) 
# Define region to plot
minregion <- max(1, peak_center - zoom_range  ) # max(1, peak_center - zoom_range  - zoom_range/2
maxregion <- peak_center + zoom_range # peak_center + zoom_range/2

region_pg <- pgParams(
  assembly = "hg38",
  chrom = chrom,
  chromstart = minregion,
  chromend = maxregion,
  x = x_start,
  width = width
)

# Load LD ------------------------------------------------------------
ld_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld0"
ld_file <- file.path(ld_dir, paste0(snp_id, ".ld"))
ld_full <- if (file.exists(ld_file)) fread(ld_file) else data.frame()

if ("SNP_B" %in% colnames(ld_full)) {
  ld_calc <- ld_full |> select(SNP_B, R2) |> rename(snp = SNP_B)
} else if ("snp" %in% colnames(ld_full)) {
  ld_calc <- ld_full |> select(snp, R2)
} else {
  ld_calc <- data.frame(snp = character(), R2 = numeric())
}

# Prepare QTL data ---------------------------------------------------
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
    rename(p = "PValue") |>
    filter(!is.na(p)) |>
    select(chrom, pos, p, snp, R2, LDgrp)|>
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

#-----------------------------------------------------------
# Layout parameters
#-----------------------------------------------------------
track_gap <- 0.125
track_height <- 0.5
manhattan_height <- 0.7
current_y <- y_start

#-----------------------------------------------------------
# PBS & FNF Manhattan plots
#-----------------------------------------------------------
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
    snpHighlights = data.frame(snp = snp_id, pch = 24, cex = 0.75, col = "black")
  )
  plotText(label = paste(cond, "caQTL"),
            x = x_start + 0.1, y = y_offset  - 0.05,
            fontsize = 7, fontfamily = "Helvetica", just = c("left", "top"))
  annoYaxis(plot = manhattan, at = seq(0, ylim_pg, 2), axisLine = TRUE, fontsize = 6)
  current_y <- current_y + manhattan_height + track_gap 
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
  
# Add legend
# legend_x <- x_start + width - 0.65
# legend_y <- y_start + height - 2.5
# plotLegend(
# legend = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
# fill = c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"),
# x = legend_x, 
# y = legend_y, 
# width = 0.1, 
# height = 0.35, 
# border = FALSE,
# fontsize = 6
# )

# Add SNP ID label
plotText(
label = snp_id,
x = x_start + width / 2, 
y = y_start - 0.1,
just = c("center", "top"),
fontfamily = "Helvetica", 
fontsize = 8
)

#-----------------------------------------------------------
# ATAC signal tracks
#-----------------------------------------------------------
if (add_atac_signal) {
  ATAC_signals <- plotMultiSignal(
    data = list(pbs_atac, fnf_atac),
    params = region_pg,
    y = current_y-0.1,
    height = track_height,
    linecolor = c(yl_gn_bu[8], yl_gn_bu[8]),
    fill = c(yl_gn_bu[8], yl_gn_bu[8]),
    default.units = "inches",
    gapdistance = 0.04
  )
  plotText(label = "ATAC", x = x_start - 0.3, y = current_y + track_height / 2 - 0.02,
            rot = 90, fontsize = 8, just = "center", fontfamily = "Helvetica", fontcolor = yl_gn_bu[8])
  plotText(label= "PBS", x = x_start-0.05, y = current_y-0.02,
        fontsize=6, fontfamily= "Helvetica", fontcolor = yl_gn_bu[8])
  plotText(label= "FN-f", x = x_start-0.05, y = current_y + track_height/2 -0.02,
        fontsize=6, fontfamily= "Helvetica", fontcolor = yl_gn_bu[8])
  current_y <- current_y + track_height
}

annoHighlight(
  plot= ATAC_signals[[1]],
  chrom = chrom, 
  chromstart = peak_start-500,
  chromend = peak_end+500,
  y = current_y-track_height-0.05 , height = 1.45,
  default.units = "inches", just = c("left", "top")
)

#-----------------------------------------------------------
# RNA signal tracks
#-----------------------------------------------------------
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

#-----------------------------------------------------------
# Gene track annotation
#-----------------------------------------------------------
gene_plot <- plotGenes(
  params = region_pg,
  y = current_y,
  height = 0.4,
  fontsize = 6
)
annoGenomeLabel(
  plot = gene_plot,
  params = region_pg,
  fontsize = 6,
  y = current_y + 0.4
)
}



# plot_caqtl_manhattan <- function(peak_id, snp_id, primary_dataset, 
#                                   x_start, y_start, width, height, 
#                                   zoom_range = 200000) {
  
# # Validate input
# if (!primary_dataset %in% c("PBS", "FNF")) {
# stop("Primary dataset must be either 'PBS' or 'FNF'")
# }

# # Select the appropriate dataset
# highConf_resQtL <- if(primary_dataset == "PBS") pbs_highconf else fnf_highconf

# # Filter data for the specified peak-SNP pair
# test_boxplotInfo <- highConf_resQtL |> 
# dplyr::filter(peak == peak_id, snp == snp_id)

# if (nrow(test_boxplotInfo) == 0) {
# stop("Peak-SNP pair not found")
# }

# # Extract chromosome and position
# peak_parts <- unlist(strsplit(peak_id, "_"))
# chrom <- peak_parts[1]
# peak_start <- as.numeric(peak_parts[2])
# peak_end <- as.numeric(peak_parts[3])

# # Load LD data
# ld_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld0"
# ld_file <- file.path(ld_dir, paste0(snp_id, ".ld"))

# # Load full chromosome LD
# ld_full <- fread(ld_file)

# if ("SNP_B" %in% colnames(ld_full)) {
#     ld_calc <- ld_full |> dplyr::select(SNP_B, R2) |> dplyr::rename(snp = SNP_B)
# } else if ("snp" %in% colnames(ld_full)) {
#     ld_calc <- ld_full |> dplyr::select(snp, R2)
# } else {
#     ld_calc <- data.frame(snp = character(), R2 = numeric())
# }

# pbs_qtl_region <- pbs_rasqual |> filter(Chromosome == chrom, peak == peak_id)
# fnf_qtl_region <- fnf_rasqual |> filter(Chromosome == chrom, peak == peak_id)

# # Join with LD
# leftjoin_pbs <- inner_join(ld_calc, pbs_qtl_region, by = "snp")
# leftjoin_fnf <- inner_join(ld_calc, fnf_qtl_region, by = "snp")
  
# # Prepare locus plot data
# prepare_locus_plot <- function(data) {
#     data |>
#         dplyr::mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) |>
#         dplyr::mutate(LDgrp = addNA(LDgrp)) |>
#         separate(snp, into = c("chrom", "pos", "ref", "alt"), sep = ":", remove = FALSE) |>
#         dplyr::mutate(pos = as.numeric(pos)) |>
#         rename(p = "PValue") |>
#         select("chrom", "pos", "p", "snp", "R2", "LDgrp") |>
#         mutate(
#         LDgrp = factor(LDgrp, 
#                         levels = c(NA, "(0.8,1]", "(0.6,0.8]", "(0.4,0.6]", 
#                                 "(0.2,0.4]", "(0,0.2]"), 
#                         ordered = TRUE)
#         ) |>
#         filter(!is.na(LDgrp)) |>
#         arrange(desc(LDgrp)) |>
#         data.frame()
# }
  
# pbs_locus_plot <- prepare_locus_plot(leftjoin_pbs)
# fnf_locus_plot <- prepare_locus_plot(leftjoin_fnf)

# # ========================================
# # Set plot region
# # ========================================
# minregion <- max(1, peak_start - zoom_range)
# maxregion <- peak_end + zoom_range

# region_pg <- pgParams(
# assembly = "hg38", 
# chrom = chrom,
# chromstart = minregion,
# chromend = maxregion,
# x = x_start, 
# width = width
# )

# # Calculate y-axis limits
# pbs_ylim <- ceiling(max(-log10(pbs_locus_plot$p), na.rm = TRUE)) + 2
# fnf_ylim <- ceiling(max(-log10(fnf_locus_plot$p), na.rm = TRUE)) + 2
# ylim_pg <- max(pbs_ylim, fnf_ylim, 5)
   
# # # Set plot region
# # minregion <- peak_start - zoom_range
# # maxregion <- peak_end + zoom_range
  
# # region_pg <- pgParams(
# #         assembly = "hg38", 
# #         chrom = chrom,
# #         chromstart = minregion,
# #         chromend = maxregion,
# #         x = x_start, 
# #         y = y_start, 
# #         width = width, 
# #         height = height
# # )
  
# # # Calculate y-axis limits
# # pbs_ylim <- ceiling(max(log10(pbs_locus_plot$p) * -1)) + 1
# # fnf_ylim <- ceiling(max(log10(fnf_locus_plot$p) * -1)) + 1
# # ylim_pg <- max(pbs_ylim, fnf_ylim)

# # ========================================
# # Track current Y position
# # ========================================
# track_gap <- 0.1
# track_height <- 0.8
# current_y <- y_start
# manhattan_height <- 1.5
# # Function to create Manhattan plot
# create_manhattan_plot <- function(data, y_offset, label) {
# plot_height <- (height - 0.2) / 2
    
# locus_plot <- plotManhattan(
#     data = data,
#     params = region_pg,
#     range = c(0, ylim_pg),
#     fill = colorby("LDgrp",
#                 palette = colorRampPalette(c("#DD3931", "#EEA741", 
#                                                 "#499A53", "#98CDED", "#262C74"))),
#     y = y_start + y_offset, 
#     height = plot_height,
#     snpHighlights = data.frame(snp = snp_id,
#                                 pch = c(24),
#                                 cex = c(0.75),
#                                 col = c("black"))
# )

    
# # Add y-axis
# annoYaxis(plot = locus_plot, at = seq(0, ylim_pg, 2),
#             axisLine = TRUE, fontsize = 6)

# # Add label
# plotText(
#     label = paste(label, "caQTL"),
#     x = x_start + 0.1, 
#     y = y_start + y_offset,
#     just = c("left", "top"),
#     fontfamily = "Helvetica", 
#     fontsize = 6
# )

# return(locus_plot)
# }
  
# # Create both PBS and FNF plots
# pbs_plot <- create_manhattan_plot(pbs_locus_plot, 0, "PBS")
# fnf_plot <- create_manhattan_plot(fnf_locus_plot, (height - 0.2) / 2 + 0.2, "FN-f")

# # Add y-axis label
# plotText(
# label = "-log10(p-value)",
# x = x_start - 0.3, 
# y = y_start + height / 2,
# rot = 90, 
# fontsize = 6, 
# just = "center",
# default.units = "inches", 
# fontfamily = "Helvetica"
# )
  
# # Add legend
# legend_x <- x_start + width - 0.5
# legend_y <- y_start + height - 2
# plotLegend(
# legend = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
# fill = c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"),
# x = legend_x, 
# y = legend_y, 
# width = 0.1, 
# height = 0.35, 
# border = FALSE,
# fontsize = 6
# )
  
# Add SNP ID label
# plotText(
# label = snp_id,
# x = x_start + width / 2, 
# y = y_start - 0.1,
# just = c("center", "top"),
# fontfamily = "Helvetica", 
# fontsize = 8
# )

# #-------------------------------------------------------
# #ATAC signal track
# #-------------------------------------------------------
# if (add_atac_signal) {
# ATAC_signals <- plotMultiSignal(
#     data = list(pbs_atac, fnf_atac),
#     params = region_pg,
#     y = current_y,
#     height = track_height,
#     linecolor = c(yl_gn_bu[3], yl_gn_bu[6]),
#     fill = c(yl_gn_bu[3], yl_gn_bu[6]),
#     default.units = "inches",
#     gapdistance = 0.02
# )
# plotText(label = "ATAC", x = x_start - 0.3, y = current_y + track_height / 2,
#             rot = 90, fontsize = 8, just = "center", fontfamily = "Helvetica")
# current_y <- current_y + track_height + track_gap
# }

# # ========================================
# # RNA Signal Tracks
# # ========================================


# # Add peak region track
# plotgenes <- plotGenes(
# params = region_pg,
# y = y_start + height + 0.1,
# height = 0.4,
# fontsize = 6
# )

# annoGenomeLabel(
# plot = plotgenes, 
# params = region_pg, 
# fontsize = 6,
# y = y_start + height + 0.55
# )
# }


