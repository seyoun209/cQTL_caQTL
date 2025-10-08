#Visualize the differentially accessible peaks from the ATAC-seq analysis
library(ggplot2)
library(ggrepel)
library(data.table)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(DESeq2)
library(httpgd)
library(ChIPseeker)
library(plotgardener)

source("/work/users/s/e/seyoun/cQTL_caQTL/scripts/atac_workflow/utils/utils_diff_atac.r")

#----- setup working directories ------
setwd("/work/users/s/e/seyoun/cQTL_caQTL/atac_output")
base_path <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/wasp"
out_data_dir <- file.path(base_path, "diff_atac/condition/data")
out_plot_dir <- file.path(base_path, "diff_atac/condition/plots")

load(file.path(out_data_dir, "02_cqnnorm_deseq2_re.RData")) # atacdds,atac_res,atac_res_Shrink
load(file.path(out_data_dir,  "03_deseq2_norm_re.RData")) # vsd, rld,normCounts
load(file.path(out_data_dir, "00_macs2_wasp_prep_dseq2data.RData")) # se, macs2_gr, counts_macs2, meta_final
load(file.path(out_data_dir,"04_atac_diff_significant_peaks.RData")) # diff_atac_sig, gained, lost, static,atac_res_Shrink_df

#----- Figure 1: Differentially Accessible Peaks -----

rownames(normCounts) <- paste0(macs2_gr$seqnames, ":", macs2_gr$start, "-", macs2_gr$end)
