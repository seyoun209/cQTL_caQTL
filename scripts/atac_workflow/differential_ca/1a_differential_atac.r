# Seyoun Byun
# date started: 09/21/2025
# date edited: 09/24/2025
# Deseq2 ATAC-seq using macs2 (regular)
#--------------------------------------------------------
#--------------------------------------------------------
source("/work/users/s/e/seyoun/cQTL_caQTL/scripts/atac_workflow/utils/utils_diff_atac.r")

#----- setup working directories ------
setwd("/work/users/s/e/seyoun/cQTL_caQTL/atac_output")
base_path <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
out_data_dir <- file.path(base_path, "diff_atac/condition/data")
out_plot_dir <- file.path(base_path, "diff_atac/condition/plots")
bamqc_dir <- file.path(base_path, "ataqv")


# Create the directory if it doesn't exist
if (!dir.exists(out_data_dir)) {
  dir.create(out_data_dir, recursive = TRUE)
}
if (!dir.exists(out_plot_dir)) {
  dir.create(out_plot_dir, recursive = TRUE)
}

#----- Loading data ------
macs2_df <- read_and_clean_peaks("peaks/counts/all_peaks_combined_counts.txt")

#Make the meta for samples
# sample config file generate
config <- yaml::read_yaml("../config/ATAC_config.yaml")
samplesheet <- fread(paste0("/work/users/s/e/seyoun/cQTL_caQTL/",config$samplesheet))
donorsheet <- fread("/work/users/s/e/seyoun/cQTL_caQTL/DonorInfo.txt", header=TRUE)

#Check
donorsheet[grep("\\(ankle\\)|\\(femur\\)", donorsheet$Donor), ]
donorsheet <- donorsheet[!grepl("\\(femur\\)", Donor)]
donorsheet[, Donor := gsub(" \\(ankle\\)", "", Donor)]

#verify
#donorsheet[grep("AM7778", donorsheet$Donor), ]

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
macs2_gr <- makeGRangesFromDataFrame(macs2_df)
macs2_gr$peakID <- macs2_df$Geneid
counts_macs2 <- as.matrix(macs2_df[,7:ncol(macs2_df)])
se <- SummarizedExperiment(assays=list(counts = counts_macs2),
                               rowRanges=macs2_gr, colData=meta_final)

save(se, macs2_gr, counts_macs2, meta_final, file = file.path(out_data_dir, "00_macs2_prep_dseq2data.RData"))


#----- QC checks -----
##  GC bias check
### Calculate peak widths
peakwidths = width(macs2_gr)

### Calculate GC content
peakseqs = getSeq(BSgenome.Hsapiens.UCSC.hg38,
                  seqnames(macs2_gr),
                  start(macs2_gr),
                  end(macs2_gr))
GCcontent = letterFrequency(peakseqs, "GC", as.prob=TRUE) |> as.numeric()

### Plot GC content distribution
gc_df <- data.frame(GCcontent = GCcontent)
mean_gc <- mean(gc_df$GCcontent, na.rm = TRUE)
median_gc <- median(gc_df$GCcontent, na.rm = TRUE)

QC_GC_histPlot <- ggplot(gc_df, aes(x = GCcontent)) +
  geom_histogram(bins = 500, fill = "grey", color = "grey") +
  geom_vline(xintercept = mean_gc, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = median_gc, color = "grey50", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = mean_gc + 0.021, y = 2300, 
   label = paste0("Mean = ", round(mean_gc, 3)),  color = "black", size = 2.5) +
  annotate("text", x = median_gc + 0.021, y = 2250, 
  label = paste0("Median = ", round(median_gc, 3)),  color = "grey50", size = 2.5) +

  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  
  labs(
    title = "Distribution of GC Content in Peaks",
    x = "GC Content",
    y = "Frequency"
  ) +
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

save(QC_GC_histPlot, file = file.path(out_plot_dir, "QC_GC_histPlot.rda"))

# Get mapped reads for samples

ataqv_files <- list.files(
  path = bamqc_dir,
  pattern = "\\.ataqv\\.txt$",
  full.names = TRUE)
ataqv_files <- ataqv_files[!grepl("Femur|replicate", ataqv_files)]
mapped_reads_df <- get_mapped_reads(ataqv_files)

# CQN normalization
cqn.counts <- cqn(counts_macs2,
                  lengths = peakwidths,
                  x = GCcontent,
                  sizeFactors = mapped_reads_df$MappedReads,
                  sqn = TRUE)

cqnOffset <- cqn.counts$glm.offset
cqnNormFactors = exp(cqnOffset)
cqnNormFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))

# save all the QC files
save(cqn.counts, cqnNormFactors,GCcontent,peakwidths,
     file = file.path(out_data_dir, "01_macs2_cqn_normFactors.RData"))

#------ DESEQ2 analysis --------------------------------------------------
#  Building deseq2 matrix

atacdds <- DESeqDataSetFromMatrix(countData = counts_macs2,
                              colData = meta_final,
                              rowRanges=macs2_gr,
                              design = ~Donor+Condition)

# Apply CQN normalization factors
normalizationFactors(atacdds) <- cqnNormFactors

keep <- rowSums(counts(atacdds) >= 10) >= ceiling(nrow(colData(se))*0.10)
atacdds  <- atacdds[keep,]

# Run DESeq2
atacdds <- DESeq(atacdds)

# Get results
atac_res <- results(atacdds)
print(summary(atac_res,alpha=0.05, na.rm=TRUE))
atac_res_Shrink <- lfcShrink(atacdds, coef='Condition_FNF_vs_CTL', type="apeglm",format =  "GRanges")
atac_res_Shrink$peakID <- rowData(atacdds)$peakID

save(atacdds,atac_res,atac_res_Shrink,
     file= file.path(out_data_dir, "02_cqnnorm_deseq2_re.RData"))

# Plotting----------------------------------------------------------
#   norm counts 

normCounts <- counts(atacdds, normalized = TRUE)
vsd <- vst(atacdds, blind=TRUE)
rld <- rlog(atacdds, blind=TRUE)

save(vsd, rld,normCounts,
 file = file.path(out_data_dir, "03_deseq2_norm_re.RData"))


# Calculate PCA on vst-transformed data
pca_vsd_re <- prcomp(t(assay(vsd)))

# Create PCA data frame
pca_df <- as.data.frame(pca_vsd_re$x)
variance_explained <- summary(pca_vsd_re)$importance[2, ]
pca_df$Donor <- colData(vsd)$Donor
pca_df$Condition <-colData(vsd)$Condition
peak_count <- nrow(assay(vsd))

PCA_vsd_scattPlot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Condition)) +
  geom_point(size = 2) +
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

# correlation plot --------------------

# Create PCA data frame
pca_df <- as.data.frame(pca_vsd_re$x)
pca_vsd_meta <- merge(pca_df, colData(vsd), by="row.names")

pbs_corr_plot <- corr_plot(pca_vsd_meta, condition = "CTL", title = "PBS-Correlations")
fnf_corr_plot <- corr_plot(pca_vsd_meta, condition = "FNF", title = "FNF-Correlations")


# MA plot------------------------------

plotMA(atac_res, ylim=c(-2,2))
atac_res_df <- atac_res |> as.data.frame()
atac_res_df$dot_class <- "NS"
atac_res_df$dot_class[atac_res_df$padj < 0.05 & atac_res_df$log2FoldChange > 0] <- "Up"
atac_res_df$dot_class[atac_res_df$padj < 0.05 & atac_res_df$log2FoldChange < 0] <- "Down"
atac_mdsPlot <- ggplot(atac_res_df, aes(x = log10(baseMean), y = log2FoldChange,color=dot_class)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c(NS = "grey60", Up = "#FFC200", Down = "#1775AE")) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
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
        y = expression("Log"[2]*" Fold change FN-f / PBS"))
#   coord_cartesian(ylim = c(-8,9))

ggsave(atac_mdsPlot, filename = file.path(out_plot_dir, "atac_mdsPlot.pdf"),
       width = 7, height = 10, units = "in", bg = "transparent")


#----------------------------------------------------------------

#----- Differential ATAC Gained, Lost, static Chondrocyte -------
## ---------- 0) Params (easy to tweak) ----------
LFC_THRESH <- 1.5
FDR_THRESH <- 0.05
ZLIM       <- 1.5      # fixed visual scale (±ZLIM)
DO_CLIP <- FALSE

## ---------- 1) Load DE results; call gained/lost/static ----------
load(file.path(out_data_dir, "02_cqnnorm_deseq2_re.RData"))
# expects: atacdds (DESeq2 object), atac_res_Shrink (DE results with shrinkage)

atac_res_Shrink <- atac_res_Shrink[!is.na(atac_res_Shrink$padj)]


condition <- colData(atacdds)$Condition
baseMean_CTL <- rowMeans(normCounts[, condition == "CTL", drop = FALSE])
baseMean_FNF <- rowMeans(normCounts[, condition == "FNF", drop = FALSE])

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


write.table(gained, file = file.path(out_data_dir, "chon_gainedFNF_specific.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lost, file = file.path(out_data_dir, "chon_lostFNF_specific.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(atac_res_Shrink_df, file = file.path(out_data_dir, "chon_atac_allPeaks.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
