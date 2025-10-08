library(data.table)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggvenn)
library(scales)
library(colorspace)
library(ggtext)
library(ggforce)
library(httpgd)
library(plotgardener)
library(colorspace)
library(mariner)
library(InteractionSet)

#------------------------------------------------------------
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

load(file.path(response_dir, "rasqual_sig_caqtl.RData")) # pbs_rasqual_sig, fnf_rasqual_sig

# Load eigenMT results
chromosomes <- 1:22
conditions <- c("pbs", "fnf")
window_kb <- 25

# eigenMT significant directories
pbs_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "pbs", "pc0")
fnf_dir <- file.path(base_dir, "caQTL/eigenMT/combined_results", paste0("window_", window_kb,"kb"), "fnf", "pc0")


# eQTL data
eqtl_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/eqtl/"
pbs_eqtl <- fread(paste0(eqtl_dir,"CTL_PEER_k20_genoPC_perm1Mb_sig_rsID_signalRanges.csv"))
fnf_eqtl <- fread(paste0(eqtl_dir,"FNF_PEER_k22_genoPC_perm1Mb_sig_rsID_signalRanges.csv"))

# sQTL data

# hic data
hic_processed_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/caQTL_processing/CQTL_AI/output/hic/00_data"
hic_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/hic_caQTL/dietJuicer/output"
load(file.path(hic_processed_dir, "diff_loop_gi.Rdata")) # diff_loop_gi, this is gi interaction clean
load(file.path(hic_processed_dir, "all_sig_loops.Rdata")) # sig_loops, this is sig loops padj < 0.1 and no log2 foldchange
load(file.path(hic_processed_dir, "fnf_gained_loops.Rdata")) # fnf_gained only sig for the fnf
load(file.path(hic_processed_dir, "pbs_gained_loops.Rdata")) # pbs_gained, fnf_lost for the fnf but I labeld wrong as the pbs gained(not gained but fnf lost

#make TXDB for the v49
#This is only one time
#txdb <- makeTxDbFromGFF(file="/proj/phanstiel_lab/Reference/human/hg38//annotations/gencode.v49.primary_assembly.annotation.gtf")
#saveDb(x=txdb, file = "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb")
txdb <- loadDb("/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb")

txdb_genes <- genes(txdb)
promoters <- promoters(txdb_genes, upstream =2500, downstream =500)

#---- sig SNPS 
pbs_caSNP <- unique(pbs_rasqual_sig$snp)
fnf_caSNP <- unique(fnf_rasqual_sig$snp)
total_caQTL_unique <- unique(c(pbs_caSNP, fnf_caSNP))

write.table(total_caQTL_unique,
    file.path(base_dir, "geno/sigcaQTLtotal.list"),
            col.names = F,
            row.names=F,
            quote=F)

uniqueIds <- data.table(V1 = unique(total_caQTL_unique))
uniqueIds[, chr := sub(":.*", "", V1)]
uniqueIds[, chr := gsub("^chr", "", chr)]

# Write out one file per chromosome
out_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/caQTL/data/sig_response_bychr"
#dir.create(out_dir, showWarnings = FALSE)

for (c in unique(uniqueIds$chr)) {
  out_file <- file.path(out_dir, paste0("chr", c, "_sig_response_snps.txt"))
  writeLines(uniqueIds[chr == c, V1], out_file)
}

# For now I will calculate LD with the 1000g but later I will do it with out samples.
# Promoter region of the caQTL peaks

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

# PBS GR:
pbs_peak_GR <- make_caQTL_GRanges(pbs_rasqual_sig$peak)
# FNF GR:
fnf_peak_GR <- make_caQTL_GRanges(fnf_rasqual_sig$peak)

# find overlapping promoters
overlap_pbscaqtlpeak_promoters <- findOverlaps(pbs_peak_GR, promoters)
overlap_fnfcaqtlpeak_promoters <- findOverlaps(fnf_peak_GR, promoters)

pbs_caPeaks_promoter <- pbs_rasqual_sig[queryHits(overlap_pbscaqtlpeak_promoters),] # 3301
fnf_caPeaks_promoter <- fnf_rasqual_sig[queryHits(overlap_fnfcaqtlpeak_promoters),] # 3583


# Get LD variants for PBS and FNF caQTLs
pbs_ld_variants <- get_caqtl_ld_variants(pbs_caPeaks_promoter, 
                                          ld_dir = "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05", 
                                          r2_threshold = 0.5)

fnf_ld_variants <- get_caqtl_ld_variants(fnf_caPeaks_promoter, 
                                          ld_dir = "/work/users/s/e/seyoun/cQTL_caQTL/atac_output/geno/ld_1000g_EUR/ld05", 
                                          r2_threshold = 0.5)


