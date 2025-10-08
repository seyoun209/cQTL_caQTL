library(plotgardener)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(yaml)
library(ggvenn)
library(scales)
library(colorspace)
library(ggtext)
library(ggforce)
library(dplyr)
library(data.table)
library(stringr)
library(httpgd)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

#--------------------------------------------------------
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

# Load eigenMT results
chromosomes <- 1:22
conditions <- c("pbs", "fnf")
window_kb <- 25

# eigenMT significant directories
pbs_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "pbs", "pc0")
fnf_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "fnf", "pc0")
load(file.path(response_dir, "rasqual_sig_caqtl.RData")) # pbs_rasqual_sig, fnf_rasqual_sig


#--------------------------------------------------------
# Make a venndiagram for the caQTLs
# Calculate unique Peak_SNP IDs in each dataset
pbs_peaks <- unique(pbs_rasqual_sig$Peak_SNP)
fnf_peaks <- unique(fnf_rasqual_sig$Peak_SNP)

# Count numbers for each set and their overlap
total_PBS <- length(pbs_peaks)
total_FNF <- length(fnf_peaks)
overlap <- length(intersect(pbs_peaks, fnf_peaks))
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap

# Calculate a scaling factor (for radius) based on the maximum count
max_count <- max(total_PBS, total_FNF)

# Build the Venn diagram plot
caQTL_venn_plot <- ggplot() +
  # PBS circle (placed at y = 0.6)
  geom_circle(aes(x0 = 0, y0 = 0.6, r = sqrt(total_PBS/max_count)),
              fill = "#BFDDFF", color = NA, alpha = 0.5) +
  # FNF circle (placed at y = -0.6)
  geom_circle(aes(x0 = 0, y0 = -0.6, r = sqrt(total_FNF/max_count)),
              fill = "#FFDDA2", color = NA, alpha = 0.5) +
  # Label for PBS-only
  geom_text(aes(x = 0, y = 0.7, label = only_PBS), size = 3) +
  # Label for FNF-only
  geom_text(aes(x = 0, y = -0.7, label = only_FNF), size = 3) +
  # Label for overlap (centered between the circles)
  geom_text(aes(x = 0, y = -0.075, label = overlap), size = 3) +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

# Display the plot
print(caQTL_venn_plot)
save(caQTL_venn_plot, file = file.path(plot_dir, "caQTL_venn_plot.rds"))

# venn for the response caQTLs
sig_response_qtl <- fread(file.path(response_dir, "significant_response_caQTLs_FDR10.tsv"))

# Compute counts
only_PBS_high <- length(pbs_highconf$Peak_SNP)
only_FNF_high <- length(fnf_highconf$Peak_SNP)
overlap_high <- length(shared_pbs$Peak_SNP == shared_fnf$Peak_SNP)
total_PBS_high <- length(unique(response_caQTL_pbs$Peak_SNP))
total_FNF_high <- length(unique(response_caQTL_fnf$Peak_SNP))

# fOR Sacling the circles, use the maximum count among the three categories
max_count_high <- max(total_PBS_high, total_FNF_high, overlap_high)

response_high_venn_plot <- ggplot() +
  # PBS circle (for high confidence, placed at y = 0.6)
  geom_circle(aes(x0 = 0, y0 = 0.6, r = sqrt(total_PBS_high/max_count_high)),
              fill = "#5fa3f1", color = NA, alpha = 0.5) +
  # FNF circle (for high confidence, placed at y = -0.6)
  geom_circle(aes(x0 = 0, y0 = -0.6, r = sqrt(total_FNF_high/max_count_high)),
              fill = "#f4b240", color = NA, alpha = 0.5) +
  # Label for PBS-only high confidence
  geom_text(aes(x = 0, y = 0.7, label = only_PBS_high), size = 3) +
  # Label for FNF-only high confidence
  geom_text(aes(x = 0, y = -0.7, label = only_FNF_high), size = 3) +
  # Label for overlap (centered between the circles)
  geom_text(aes(x = 0, y = -0.075, label = overlap_high), size = 3) +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

save(response_high_venn_plot, file = file.path(plot_dir, "response_highconf_caQTL_venn_plot.rds"))

#--------------------------------------------------------
# Chromhmm (e049) barplot

#function
calc_chromhmm_enrichment <- function(peak_vec, chromhmm_filepath, states,
                                       bsgenome = BSgenome.Hsapiens.UCSC.hg38) {

  # Prepare peaks GRanges from peak_vec (assumed format: "chr_peakstart_peakend_...")
  peaks <- unique(peak_vec)
  peaks <- na.omit(peaks)
  peak_split <- str_split_fixed(peaks, "_", 4)
  chr <- peak_split[,1]
  start <- as.numeric(peak_split[,2])
  end <- as.numeric(peak_split[,3])

  peak_df <- data.frame(CHR = chr, START = start, END = end, stringsAsFactors = FALSE)
  peak_GR <- makeGRangesFromDataFrame(peak_df)

  # Calculate genome length from the provided BSgenome
  genomelength <- sum(as.numeric(seqlengths(bsgenome)), na.rm = TRUE)

  # Load ChromHMM file and set appropriate column names
  chromhmm_temp <- fread(chromhmm_filepath, skip = 1, header = FALSE)
  colnames(chromhmm_temp) <- c("chr","start","end","name","score","strand","ThickStart","ThickEnd","ItemRgb")  

  # Initialize output data frame for enrichment results
  Output <- as.data.frame(matrix(NA, nrow = length(states), ncol = 6))
  rownames(Output) <- states
  colnames(Output) <- c("caQTL_Enrichment", "PercentOverlaps", "nlog10pval",
                        "Enrichment_Score", "Overlap_Count", "TotalPeaks")
  # Loop over each ChromHMM state
  for (j in seq_along(states)) {
    # Subset ChromHMM data to current state
    thisstate <- chromhmm_temp[chromhmm_temp$name == states[j], ]

    # Create GRanges for the state regions  
    thisstateGR <- GRanges(seqnames = thisstate$chr,
                           ranges = IRanges(start = thisstate$start, end = thisstate$end))
    # Filter to chromosomes present in the peak data
    thisstateGR <- thisstateGR[as.character(seqnames(thisstateGR)) %in% as.character(seqlevels(peak_GR))]

    # Collapse overlapping regions if any
    TSGR <- reduce(thisstateGR, ignore.strand = TRUE)

    # Find overlaps between peaks and state regions
    hits <- findOverlaps(peak_GR, thisstateGR)
    s_binom <- length(unique(queryHits(hits)))  # overlapping peak count
    n_binom <- length(peak_GR)                   # total peaks

    # Fraction of the genome covered by this state
    p_binom <- sum(width(TSGR)) / genomelength

    # Binomial test: probability of observing >= s_binom overlaps
    pval <- pbinom(s_binom - 1, n_binom, p_binom, lower.tail = FALSE, log.p = FALSE)
    nlog10pval <- -log10(pval)

    # Enrichment score calculation:
    # (observed bp overlap / genome length) divided by
    # [ (total bp in peaks/genome length) * (state fraction) ]
    numbpsoverlap <- sum(width(intersect(peak_GR, TSGR)))
    expected <- (sum(width(peak_GR)) / genomelength) * p_binom
    enrichment_score <- (numbpsoverlap / genomelength) / expected

    # Store the calculated metrics
    Output[j, "caQTL_Enrichment"] <- pval
    Output[j, "PercentOverlaps"] <- s_binom / n_binom
    Output[j, "nlog10pval"] <- nlog10pval
    Output[j, "Enrichment_Score"] <- enrichment_score
    Output[j, "Overlap_Count"] <- s_binom
    Output[j, "TotalPeaks"] <- n_binom
  }

  return(Output)
}

## Define ChromHMM states
states <- c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk",
            "6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv",
            "11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies")

pbs_Output <- calc_chromhmm_enrichment(pbs_rasqual_sig$peak,
                                       "/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/tools/chromhmm/E049_15_coreMarks_hg38lift_dense.bed.gz",
                                       states)
pbs_Output$fdr <- p.adjust(pbs_Output$caQTL_Enrichment, method = "fdr")

fnf_Output <- calc_chromhmm_enrichment(fnf_rasqual_sig$peak,
                                       "/proj/phanstiel_lab/Data/processed/CQTL/sqtl/CQTL_sQTL/tools/chromhmm/E049_15_coreMarks_hg38lift_dense.bed.gz",
                                       states)
fnf_Output$fdr <- p.adjust(fnf_Output$caQTL_Enrichment, method = "fdr")

pbs_Output$Condition <- "PBS"
pbs_Output$State <- rownames(pbs_Output)
fnf_Output$Condition <- "FNF"
fnf_Output$State <- rownames(fnf_Output)
# Combine
combined <- rbind(pbs_Output, fnf_Output)

# Make barplot

combined$Group <- dplyr::case_when(
  combined$State %in% c("1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk") ~ "Promoter",
  combined$State %in% c("3_TxFlnk", "4_Tx", "5_TxWk") ~ "Transcribed",
  combined$State %in% c("6_EnhG", "7_Enh", "12_EnhBiv") ~ "Enhancer",
  combined$State %in% c("13_ReprPC", "14_ReprPCWk") ~ "Polycomb",
  combined$State %in% c("9_Het") ~ "Heterochromatin",
  combined$State %in% c("8_ZNF/Rpts") ~ "ZNF/Repeats",
  combined$State %in% c("15_Quies") ~ "Quiescent",
  TRUE ~ "Other"
)

# Set the desired color mapping for the conditions
combined$Condition <- factor(combined$Condition, levels = c("PBS","FNF"))
# Create a custom label column for the conditions:
combined$ConditionLabel <- ifelse(combined$Condition == "FNF", "FN-f", "PBS")

# Define the desired order for the broader groups
group_order <- c("Promoter", "Transcribed", "Enhancer", "Heterochromatin", "ZNF/Repeats", "Polycomb", "Quiescent")
combined$Group <- factor(combined$Group, levels = group_order)


# Also ensure the State factor is ordered correctly (as defined before)
state_order <- c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk",
                 "6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv",
                 "11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies")
combined$State <- factor(combined$State, levels = state_order)

save(combined, pbs_Output, fnf_Output,file = file.path(plot_data_dir, "caQTL_chromhmm_enrichment.Rdata"))

# sig definie label based on p-values
combined <- combined %>%
  mutate(label = dplyr::case_when(
    fdr > 0.05 ~ "",
    fdr > 0.01 ~ "*",
    fdr > 0.001 ~ "**",
    TRUE ~ "***"
  ))


chromhmm_barplot <- ggplot(combined, aes(x = State, y = Enrichment_Score, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = label, y = Enrichment_Score + 0.5),
            position = position_dodge(width = 0.9), size = 3) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1, strip.position = "bottom") +
  labs(y = "Enrichment Score",
       title = "ChromHMM Enrichments for PBS vs FNF caQTL Peaks") +
  scale_fill_manual(values = c("PBS" = "#2057A7", "FNF" = "#F2BC40"),
                    labels = c("PBS", "FN-f")) +
  theme_pubr(x.text.angle = 60) +
  theme(plot.title = element_text(hjust = 0.5, size = 6),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        axis.title.y = element_text(size = 6),
        title = element_text(size=6),
        strip.placement = "outside",
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        strip.background = element_blank(),
        strip.text.x.bottom = element_text(size = 6),
        legend.key.size = unit(0.1, "cm"),
        legend.title = element_text(size = 4),
        legend.text  = element_text(size = 4))

save(chromhmm_barplot, file = file.path(plot_dir, "chromhmm_barplot.rds"))

# Combined 1st section of the caQTL figure----------------------------------------------------------
pdf(file.path(plot_dir, "caQTL_chromhmm.pdf"), width=5, height=7, bg="transparent")
pageCreate(width = 5, height = 7, showGuides = FALSE)
plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file.path(plot_dir, "caQTL_venn_plot.rds"))

plotText("Significant genetic and condition interaction effect", 
                just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8,x = 2, y = 0.35)

plotGG(caQTL_venn_plot, x = 0.63, y = 0.5, width = 1.75, height = 2.5)

plotText("PBS \nre-caQTL",
         just = c("center","top"), fontfamily = "Helvetica",fontsize = 8,x = 0.5, y = 1.3, fontcolor = "#2057A7",fontface = "bold" )

plotText("FN-f \nre-caQTL",
         just =c("center","top"), fontfamily = "Helvetica",fontsize = 8,x = 0.5, y = 2.25, fontcolor ="#F2BC40",fontface = "bold" )

plotSegments(
  x0 = c(1.5,2.5), y0 = c(0.7,0.55), x1 = 1.5, y1 = 0.55,
  default.units = "inches",
  lwd = 1, lty = 2
)

plotSegments(
  x0 = c(1.5,2.5), y0 = c(2.75,2.9), x1 = 1.5, y1 = 2.9,
  default.units = "inches",
  lwd = 1, lty = 2
)

plotGG(response_high_venn_plot, x = 2.5, y =1.25, width = 1.2, height = 1.2)

plotText("high confidence \nPBS-specific",
         just = c("center","top"), fontfamily = "Helvetica",fontsize = 7,x = 3, y = 1, fontcolor = "#2057A7",fontface="bold")

plotText("nhigh confidence \nFN-f-specific",
         just =c("center","top"), fontfamily = "Helvetica",fontsize = 7,x = 3, y = 2.5, fontcolor ="#F2BC40",fontface="bold" )

# Chromhmm plot
plotText("b", x = 0.1, y = 3.5, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
#chromhmm_barplot,
load(file.path(plot_dir, "chromhmm_barplot.rds"))
plotGG(chromhmm_barplot, x = 0.3, y = 4.0, width = 4.5, height = 3)

dev.off()
