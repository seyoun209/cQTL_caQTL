# Seyoun Byun
# date started: 08.27.2025
# date edited: 08.29.2025
# This is script to use the wasp for deseq2 ATAC-seq using macs2
#--------------------------------------------------------
#----- Libraries------
library(DESeq2)
library(Rsubread,lib.loc = "/nas/longleaf/home/seyoun/R/x86_64-pc-linux-gnu-library/4.2")
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)
library(base)
library(data.table)
library(tidyverse)
library(dplyr)
library(tibble)
library(cqn)
library(rasqualTools)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(msigdbr)
library(plotgardener)
library(tximeta)
library(readr)
library(rrvgo)
library(gridtext)
library(ggtext)
library(colorspace)
library(httpgd)

#----- setup working directories ------
setwd("/proj/phanstiel_lab/Data/processed/CQTL/caQTL_processing/CQTL_AI/output")
base_path <- "/proj/phanstiel_lab/Data/processed/CQTL/caQTL_processing/CQTL_AI"
out_data_dir <- file.path(base_path, "output/diff_wasp_deseq2/condition/data")
out_plot_dir <- file.path(base_path, "output/diff_wasp_deseq2/condition/plots")
#----- Functions ------
## Read and clean data
read_and_clean_peaks <- function(file) {
  df <- read.table(file, header=TRUE, skip=1, sep="\t")
  colnames(df) <- gsub("output.wasp.blk_filter.", "", colnames(df))
  colnames(df) <- gsub(".sorted_final.bam", "", colnames(df))
  return(df)
}

# Create the directory if it doesn't exist
if (!dir.exists(out_data_dir)) {
  dir.create(out_data_dir, recursive = TRUE)
}

if (!dir.exists(out_plot_dir)) {
  dir.create(out_plot_dir, recursive = TRUE)
}

#z-score calculation
cal_z_score <- function(x) {
  s <- sd(x)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x)) / s
}
#----- Loading data ------
macs2  <- read_and_clean_peaks("wasp_peaks/merged/allsamples_macs2_merged_counts.txt")

#Make the meta for samples
# sample config file generate (Metasamplesheet)
config <- yaml::read_yaml("../config/ATACconfig.yaml")
samplesheet <- fread(paste0("/work/users/s/e/seyoun/CQTL_AI/",config$samplesheet))
donorsheet <- fread("/work/users/s/e/seyoun/CQTL_AI/Chon_DonorInfo.txt") 

#Check
donorsheet[grep("\\(ankle\\)|\\(femur\\)", donorsheet$Donor), ]
donorsheet <- donorsheet[!grepl("\\(femur\\)", Donor)]
donorsheet[, Donor := gsub(" \\(ankle\\)", "", Donor)]

#verify
donorsheet[grep("AM7778", donorsheet$Donor), ]

merged_meta <- left_join(samplesheet, donorsheet, by = "Donor") |> 
  filter(Tissue == "Ankle") |>
  dplyr::select("Proj","Donor","Condition","Tissue","Protocol_notes","Time","Tech_Rep","Seq_Rep","TreatmentDate",
                "ATACProtocolDate","ATACLibrarySubmissionDate","Sex","Age","Race",
                "OAGradeAvg","CauseOfDeath","FragmentBatch")



#Take the replicate out
pairs_count <- merged_meta[, .N, by=.(Donor, Condition)]
duplicated_pairs <- pairs_count[N > 1, .(Donor, Condition)]

meta_data_rep_cleaned <- merged_meta[!(Donor %in% unique(duplicated_pairs$Donor)) | 
                                       (Donor %in% unique(duplicated_pairs$Donor) & Tech_Rep == 1)]

# Range of the age group
meta_data_rep_cleaned$AgeGroup <- cut(meta_data_rep_cleaned$Age,
                                      breaks=c(24, 44, 64, 84),
                                      labels=c("25-44", "45-64", "65-84"))



meta_final <- meta_data_rep_cleaned |> 
  mutate(sampleID = paste(Proj, Donor, Condition, Tissue, Protocol_notes, sep="_"))

meta_final$Donor <- factor(meta_final$Donor)
meta_final$Condition <- factor(meta_final$Condition, levels = c("CTL", "FNF"))

#Change to macs2_gr
macs2_gr <- makeGRangesFromDataFrame(macs2)
macs2_gr$peakID <- macs2$Geneid
counts_macs2 <- as.matrix(macs2[,7:ncol(macs2)])
se <- SummarizedExperiment(assays=list(counts = counts_macs2),
                               rowRanges=macs2_gr, colData=meta_final)

save(se, macs2_gr, counts_macs2, meta_final, file = file.path(out_data_dir, "chon_macs2_wasp_se.RData"))


#----- QC checks -----
## 1. GC bias check
### Calculate peak widths
peakwidths = width(macs2_gr)

### Calculate GC content
message('Calculating GC content...')
peakseqs = getSeq(BSgenome.Hsapiens.UCSC.hg38,
                  seqnames(macs2_gr),
                  start(macs2_gr),
                  end(macs2_gr))
GCcontent = letterFrequency(peakseqs, "GC", as.prob=TRUE) |> as.numeric()

### Plot GC content distribution
gc_df <- data.frame(GCcontent = GCcontent)
mean_gc <- mean(gc_df$GCcontent, na.rm = TRUE)

gc_dist_plotHist <- ggplot(gc_df, aes(x = GCcontent)) +
  geom_histogram(bins = 500, fill = "grey", color = "grey") +
  geom_vline(xintercept = mean_gc, color = "black", linetype = "dashed", size = 0.5) +
  annotate("text", x = mean_gc, y = 3200, label = paste0("Mean = ", round(mean_gc, 3)),  color = "black", size = 2.5) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Distribution of GC Content in Peaks",
    x = "GC Content",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(color = "black", size = 8),
    axis.line.x = element_line(color = "black", linewidth = 0.25),
    axis.line.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_text(color = "black", size = 8),
    plot.title = element_text(color = "black", size = 8),
    plot.subtitle = element_blank()
  ) 

save(GCcontent,gc_dist_plotHist,
     file = file.path(out_data_dir, "GCContent.RData"))


## Get library size information
message('Reading library size metrics...')
### Shell to run
### ml multiqc/1.25.1

flagstat_data <- fread("wasp/blk_filter/qc_cqtl/multiqc_data/multiqc_samtools_flagstat.txt",
                      header = TRUE,
                      stringsAsFactors = FALSE)

mapped_reads <- data.frame(
  sampleID = gsub("_stats$", "", flagstat_data$Sample),
  MappedReads = flagstat_data$mapped_passed
) %>%
  filter(sampleID %in% colnames(counts_macs2)) %>%
  dplyr::slice(match(colnames(counts_macs2), sampleID))

## Apply CQN normalization
message('Performing CQN normalization...')
cqn.counts <- cqn(counts_macs2,
                  lengths = peakwidths,
                  x = GCcontent,
                  sizeFactors = mapped_reads$MappedReads,
                  sqn = TRUE)

### Calculate normalization factors for RASQUAL
message('Calculating RASQUAL normalization factors...')
cqnOffset <- cqn.counts$offset
cqnNormFactors <- exp(-cqnOffset)  # Convert offsets to factors (note negative sign)

### Center normalization factors around 1.0 (optional but recommended)
cqnNormFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))

### Create normalized counts (adding this missing piece)
message('Calculating normalized expression values...')
normalized_counts <- cqn.counts$y + cqn.counts$offset


### Save the CQN results for inspection/QC
save(cqn.counts, GCcontent, peakwidths,cqnNormFactors,
     file = file.path(out_data_dir, "CQN_results.RData"))

#--- visualization
### Before correction

# Save smoothScatter plots as R objects (no text annotation)
pre_lm <- lm(rowMeans(log2(counts_macs2 + 1)) ~ GCcontent)
pre_slope <- coef(pre_lm)[2]


pre_scatter_plot <- function() {
  smoothScatter(GCcontent, rowMeans(log2(counts_macs2 + 1)),
        main="Before GC Correction",
        xlab="GC Content",
        ylab="Mean log2(Counts + 1)")
  abline(pre_lm, col="red", lwd=2)
  text(0.2, max(rowMeans(log2(counts_macs2 + 1)))*0.9, 
    paste("Slope =", round(pre_slope, 2)), col="red")
}

post_lm <- lm(rowMeans(normalized_counts) ~ GCcontent)
post_slope <- coef(post_lm)[2]


post_scatter_plot <- function() {
  smoothScatter(GCcontent, rowMeans(normalized_counts),
        main="After CQN Correction",
        xlab="GC Content",
        ylab="Mean Normalized Expression")
  abline(post_lm, col="red", lwd=2)
  text(0.2, max(rowMeans(normalized_counts))*0.9, 
    paste("Slope =", round(post_slope, 2)), col="red")
}

save(pre_scatter_plot, post_scatter_plot, pre_slope, post_slope, file = file.path(out_plot_dir, "GC_bias_scatterplots.RData"))



#---------------------------------------------------------------
#  Building deseq2 matrix

atacdds <- DESeqDataSetFromMatrix(countData = counts_macs2,
                              colData = meta_final,
                              rowRanges=macs2_gr,
                              design = ~Donor+Condition)

# Apply CQN normalization factors
normalizationFactors(atacdds) <- cqnNormFactors

# Run DESeq2
atacdds <- DESeq(atacdds)

# Get results
atac_res <- results(atacdds)
print(summary(atac_res,alpha=0.05, na.rm=TRUE))
atac_res_Shrink <- lfcShrink(atacdds, coef='Condition_FNF_vs_CTL', type="apeglm",format =  "GRanges")
atac_res_Shrink$peakID <- rowData(atacdds)$peakID

save(atacdds,atac_res,atac_res_Shrink,
     file= file.path(out_data_dir, "wasp_atac_cqnnorm_deseq2_re.RData"))


#----------------------------------------------------------------
# 3. PCA plot

message("Calculating PCA for RASQUAL covariates...")
# Calculate variance-stabilized transformation
normCounts <- counts(atacdds, normalized = TRUE)
vsd <- vst(atacdds, blind=TRUE)
rld <- rlog(atacdds, blind=TRUE)

save(vsd, rld,normCounts,
 file = file.path(out_data_dir, "wasp_macs2_deseq_cqn_normalized.RData"))

# Calculate PCA on vst-transformed data
pca_result <- prcomp(t(assay(vsd)))

# Create PCA data frame
pca_df <- as.data.frame(pca_result$x)
variance_explained <- summary(pca_result)$importance[2, ]
pca_df$Donor <- colData(vsd)$Donor
pca_df$Condition <-colData(vsd)$Condition


# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_var <- cumsum(var_explained)

# Plot variance explained
pdf(file.path(out_plot_dir,"varianceplot.pdf"), width=10, height=5)
par(mfrow=c(1,2))
barplot(var_explained[1:20]*100, 
        main="Variance Explained by PCs",
        xlab="Principal Component", 
        ylab="Percent Variance Explained",
        names.arg=1:20)
plot(cumulative_var[1:20]*100, 
     type="b", pch=19,
     xlab="Principal Component",
     ylab="Cumulative Percent Variance Explained",
     main="Cumulative Variance")
abline(h=75, col="red", lty=2)
dev.off()


peak_count <- nrow(assay(vsd))
pca_plot_vsd <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("CTL" = "#1775AE", "FNF" = "#FFC200")) +
  labs(
    title = paste0("Total peaks (n=", peak_count, ") with all samples"),
    x = paste0("PC1 Variance: ", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))


print(pca_plot_vsd)
# Save PCA plot
save(pca_plot_vsd,
file = file.path(out_plot_dir,"PCA_plot_vsd.rds"))


# Create PCA plots with different grouping variables using DESeq2's plotPCA
plotPCA(rld, intgroup = c("Condition", "AgeGroup"), ntop=400) + 
  ggplot2::theme(aspect.ratio = 1)
plotPCA(vsd, intgroup = c("Condition", "Sex"), ntop=400) + 
  ggplot2::theme(aspect.ratio = 1)

# Alternative PCA using rlog transformation for comparison
mat_atac <- assay(rld)
pca_comp <- prcomp(t(mat_atac))
variance_explained_rld <- summary(pca_comp)$importance[2, ]

# Create PCA data frame
pca_df_rld <- as.data.frame(pca_comp$x)
pca_df_rld$Donor <- colData(rld)$Donor
pca_df_rld$Condition <- colData(rld)$Condition

# Create PCA plot with rlog data
pca_plot_rld <- ggplot(pca_df_rld, aes(x = PC1, y = PC2, label = Donor, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("CTL" = "#1775AE", "FNF" = "#FFC200")) +
  labs(
    title = paste0("Total peaks (n=", nrow(mat_atac), ") with all samples (rlog)"),
    x = paste0("PC1 Variance: ", round(variance_explained_rld[1] * 100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained_rld[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))

print(pca_plot_rld)
save(pca_plot_rld,
     file = file.path(out_plot_dir, "PCA_plot_rld.rds"))


#----- MA plot -----

atac_maplot <- ggplot(atac_res, aes(x = log10(baseMean), y = log2FoldChange)) +
   geom_point(data = subset(atac_res, padj > 0.05), 
              color = "grey60", alpha = 0.5 )+
   geom_point(data = subset(atac_res, padj < 0.05 & log2FoldChange > 0), 
              color = "#FFC200", alpha = 0.5) +
   geom_point(data = subset(atac_res, padj < 0.05 &  log2FoldChange < 0), 
              color = "#1775AE", alpha = 0.5) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
   theme_classic() +
   theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
         plot.background = element_rect(fill = 'transparent', color = "transparent"),
         text = element_text(family = "Helvetica"),
         legend.position = "None",
         axis.text.y = element_text(color = "black", size = 6),
         axis.title.y = element_text(size = 6),  # Changed from element_markdown
         axis.title.x = element_text(size = 6),  # Changed from element_markdown
         axis.text.x = element_text(color = "black", size = 6),
         strip.background = element_blank(),
         axis.ticks = element_blank(),
         axis.line.x = element_line(linewidth = 0.25),
         strip.text = element_text(size = 8, color = "black"),
         panel.spacing = unit(0, "mm"), 
         plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
   labs(x = expression("Mean of normalized counts (log"[10]*")"),
        y = expression("Log"[2]*" Fold change FN-f / PBS"))+
   coord_cartesian(ylim = c(-8,9))


save(atac_maplot,
     file = file.path(out_plot_dir, "MA_plot_rld.rds"))


#----------------------------------------------------------------

#----- Differential ATAC Gained, Lost, static Chondrocyte -------
## ---------- 0) Params (easy to tweak) ----------
LFC_THRESH <- 1.5
FDR_THRESH <- 0.05
ZLIM       <- 1.5      # fixed visual scale (±ZLIM)
DO_CLIP <- FALSE

## ---------- 1) Load DE results; call gained/lost/static ----------
load(file.path(out_data_dir, "wasp_atac_cqnnorm_deseq2_re.RData"))
# expects: atacdds (DESeq2 object), atac_res_Shrink (DE results with shrinkage)

atac_res_Shrink <- atac_res_Shrink[!is.na(atac_res_Shrink$padj)]
normalized_counts <- counts(atacdds, normalized = TRUE)

condition <- colData(atacdds)$Condition
baseMean_CTL <- rowMeans(normalized_counts[, condition == "CTL", drop = FALSE])
baseMean_FNF <- rowMeans(normalized_counts[, condition == "FNF", drop = FALSE])

atac_res_Shrink_df <- as.data.frame(atac_res_Shrink)

diff_atac_sig <- atac_res_Shrink_df %>%
  filter(abs(log2FoldChange) > LFC_THRESH, padj < FDR_THRESH)

gained <- diff_atac_sig %>% filter(log2FoldChange > 0)
lost   <- diff_atac_sig %>% filter(log2FoldChange < 0)
static <- atac_res_Shrink_df %>%
  filter(!(peakID %in% c(gained$peakID, lost$peakID)))

atac_res_Shrink_df$class <- "static"
atac_res_Shrink_df$class[atac_res_Shrink_df$peakID %in% gained$peakID] <- "gained"
atac_res_Shrink_df$class[atac_res_Shrink_df$peakID %in% lost$peakID]   <- "lost"


write.table(gained, file = file.path(out_data_dir, "chondrocyte_FNF_specific.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lost, file = file.path(out_data_dir, "chondrocyte_PBS_specific.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(atac_res_Shrink_df, file = file.path(out_data_dir, "chondrocyte_atac_allPeaks.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)



#--------------------------------------------------------------
# Heatmap
load(file = file.path(out_data_dir, "chon_macs2_wasp_se.RData"))
load(file = file.path(out_data_dir, "wasp_macs2_deseq_cqn_normalized.RData"))

# Filter out the significant counts from the normalized counts
rownames(normCounts) <- rowData(se)$peakID
range_info <- data.frame(
  seqnames = as.character(seqnames(rowRanges(se))),
  start    = start(rowRanges(se)),
  end      = end(rowRanges(se)),
  row.names = rowData(se)$peakID
)

# keep only significant peaks in heatmap
keep_ids <- intersect(rownames(normCounts), diff_atac_sig$peakID)
norm_counts_sig <- normCounts[keep_ids, , drop = FALSE]
zscore_matrix <- t(apply(norm_counts_sig, 1, cal_z_score))

#drop rows with any NA (shouldn't happen, but be safe)
zscore_matrix <- zscore_matrix[rowSums(is.na(zscore_matrix)) == 0, , drop = FALSE]

# only the samples in the matrix, preserve their current order
sample_meta <- meta_final %>%
  filter(sampleID %in% colnames(zscore_matrix)) %>%
  select(sampleID, Condition, Sex, AgeGroup)
rownames(sample_meta) <- sample_meta$sampleID

# Relevel factors for legends (these are display labels; not used in modeling)
sample_meta$Condition <- factor(sample_meta$Condition, levels = c("CTL","FNF"))
sample_meta$Sex       <- factor(sample_meta$Sex, levels = c("M","F"))
sample_meta$AgeGroup  <- factor(sample_meta$AgeGroup, levels = c("25-44","45-64","65-84"))

# make sure sample_meta is rownamed by sampleID
if (is.null(rownames(sample_meta)) || any(rownames(sample_meta) != sample_meta$sampleID)) {
  rownames(sample_meta) <- sample_meta$sampleID
}

# named vector of conditions (bulletproof vs data.table)
cond_vec <- setNames(as.character(sample_meta$Condition), rownames(sample_meta))

ID <- colnames(zmat)
column_order <- c(ID[cond_vec[ID] == "CTL"], ID[cond_vec[ID] == "FNF"])

meta_ctl_fnf <- sample_meta |>
  mutate(ord = match(sampleID, ID),
         Condition = factor(Condition, levels = c("CTL","FNF"))) |>
  arrange(Condition, ord)
column_order <- meta_ctl_fnf$sampleID

# Annotation df via named vectors (no [i,j] on tables)
age_vec <- setNames(as.character(sample_meta$AgeGroup), rownames(sample_meta))
sex_vec <- setNames(as.character(sample_meta$Sex),       rownames(sample_meta))

ann_df <- data.frame(
  AgeGroup  = factor(age_vec[column_order], levels = c("25-44","45-64","65-84")),
  Sex       = factor(sex_vec[column_order], levels = c("M","F")),
  Condition = factor(cond_vec[column_order], levels = c("CTL","FNF")),
  row.names = column_order, stringsAsFactors = FALSE
)


## ---------- DAR row metadata & split ----------
dar_df <- atac_res_Shrink_df %>%
  select(peakID, log2FoldChange, padj) %>%
  filter(peakID %in% rownames(zscore_matrix)) %>%
  column_to_rownames("peakID")

dir_vec <- factor(ifelse(dar_df$log2FoldChange >= 0, "Gained", "Lost"),
                  levels = c("Gained","Lost"))

## ---------- Colors & (optional) clipping ----------
heat_cols <- circlize::colorRamp2(c(-ZLIM, 0, ZLIM),
                                  c("#097EA4", "black", "#BFA527"))
zmat <- if (DO_CLIP) pmax(pmin(zscore_matrix, ZLIM), -ZLIM) else zscore_matrix

## ---------- Annotations ----------
my_colour <- list(
  Condition = c("CTL" = "#2057A7", "FNF" = "#F2BC40"),
  Sex       = c("M" = "#89d1dcff", "F" = "#ef8871ff"),
  AgeGroup  = c("25-44" = "#b5d1ae", "45-64" = "#80ae9a", "65-84" = "#568b87")
)

colAnn <- HeatmapAnnotation(
  df = ann_df,
  col = my_colour,
  gap = unit(0.5, "mm"),
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    AgeGroup = list(title="AgeGroup", title_position="leftcenter",
                    title_gp=gpar(fontsize=7), labels_gp=gpar(fontsize=6),
                    grid_width=unit(2,"mm"), grid_height=unit(1,"mm"),
                    at=c("25-44","45-64","65-84"),
                    labels=c("25-44","45-64","65-84"), ncol=1),
    Condition = list(title="Condition", title_position="leftcenter",
                     title_gp=gpar(fontsize=7), labels_gp=gpar(fontsize=6),
                     grid_width=unit(2,"mm"), grid_height=unit(1,"mm"),
                     at=c("CTL","FNF"), labels=c("PBS","FNF"), ncol=1),
    Sex = list(title="Sex", title_position="leftcenter",
               title_gp=gpar(fontsize=7), labels_gp=gpar(fontsize=6),
               grid_width=unit(2,"mm"), grid_height=unit(1,"mm"),
               at=c("M","F"), labels=c("Male","Female"), ncol=1)
  )
)

## ---------- Heatmap ----------
hmap <- Heatmap(
  zmat[, column_order, drop = FALSE],
  col = heat_cols,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = column_order,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE,
  top_annotation = colAnn,
  use_raster = FALSE  # set TRUE if you want smaller files
)

## ---------- Draw & capture grobs ----------
heatmapLegend <- Legend(
  at = c(-ZLIM, 0, ZLIM),
  col_fun = heat_cols,
  border = NA,
  title_gp = gpar(fontsize = 0),
  labels_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
  legend_width = unit(4.325, "in"),
  grid_height = unit(0.11, "in"),
  direction = "horizontal"
)

heatmapGrob <- grid.grabExpr(draw(hmap, show_annotation_legend=FALSE,
                                        show_heatmap_legend=FALSE, background="transparent"))
heatmapLegendGrob <- grid.grabExpr(draw(heatmapLegend))

save(heatmapGrob,       file = file.path(out_plot_dir, "heatmapGrob.rda"))
save(heatmapLegendGrob, file = file.path(out_plot_dir, "heatmapLegendGrob.rda"))
save(zscore_matrix,     file = file.path(out_plot_dir, "zscore_matrix_raw.rda"))



#------------------------------------------------
#----- Create boxplot for the promotor regions  ------
txdb <- loadDb("/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)
promoters <- promoters(txdb_genes, upstream =1500, downstream =500)

load(file.path(out_data_dir, "wasp_atac_cqnnorm_deseq2_re.RData"))
# expects: atacdds (DESeq2 object), atac_res_Shrink (DE results with shrinkage)

atac_res_Shrink <- atac_res_Shrink[!is.na(atac_res_Shrink$padj)]
diff_atac_sig_gr <- atac_res_Shrink |> 
  filter(abs(log2FoldChange) > LFC_THRESH, padj < FDR_THRESH)

gained_gr <- diff_atac_sig_gr |> filter(log2FoldChange > 0)
lost_gr   <- diff_atac_sig_gr |> filter(log2FoldChange < 0)
static_gr <- diff_atac_sig_gr |>
  filter(!(peakID %in% c(gained_gr$peakID, lost_gr$peakID)))

atac_res_Shrink$class <- "static"
atac_res_Shrink$class[atac_res_Shrink$peakID %in% gained_gr$peakID] <- "gained"
atac_res_Shrink$class[atac_res_Shrink$peakID %in% lost_gr$peakID]   <- "lost"


# find overlapping promoters
overlap_gainedatac_promoters <- findOverlaps(gained_gr, promoters)
overlap_staticatac_promoters <- findOverlaps(static_gr, promoters)
overlap_lostatac_promoters   <- findOverlaps(lost_gr, promoters)

# Get the overlaping gene those

promoters_at_gained_atac <- promoters[subjectHits(overlap_gainedatac_promoters)]
promoters_at_static_atac <- promoters[subjectHits(overlap_staticatac_promoters)]
promoters_at_lost_atac   <- promoters[subjectHits(overlap_lostatac_promoters)]


#------------------------------------------------
# Get the genes from the RNA-seq differential expression results


#------------------------------------------------
#----- Create the motif analysis using Homer ------

#------------------------------------------------
#----- Create Fig1-----

pdf(file=file.path(out_plot_dir, "fig1_wasp_atac.pdf"),width=13,height=6,bg="transparent")  # Set background to transparent
pageCreate(width = 13, height = 6, showGuides = TRUE) 
plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file.path(out_plot_dir, "heatmapGrob.rda"))  # Load the heatmap grob
load(file.path(out_plot_dir, "heatmapLegendGrob.rda"))  # Load the heatmap legend grob

plotGG(plot = heatmapGrob, x = 0.45, y = 0.45, height = 4.5, width = 5)
plotGG(plot = heatmapLegendGrob, x = 0.55, y = 4.95,
       width = 4.325, height = 0.11)


dev.off()
