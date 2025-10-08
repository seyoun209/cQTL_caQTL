# Size factors
# Date: Sept 2025
#------------------------------------------------------
# Load libraries
library(dplyr)
library(stringr)
library(GenomicRanges)
library(rasqualTools)

#------------------------------------------------------
# Load data
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
 
expression_out_dir <- file.path(base_dir, "caQTL","rasqual_input","expression")
sizefactor_out_dir <- file.path(base_dir, "caQTL","rasqual_input","sizefactor")
peak_info_dir <- file.path(base_dir, "caQTL","rasqual_input","peakinfo") 

# se, macs2_gr, counts_macs2, meta_final
load(file.path(base_dir,"wasp","diff_atac","condition","data","00_macs2_wasp_prep_dseq2data.RData"))
# cqn.counts, cqnNormFactors,GCcontent,peakwidths
load(file.path(base_dir,"wasp","diff_atac","condition","data", "01_macs2_wasp_cqn_normFactors.RData"))

chromosomes <- paste0("chr", c(1:22))

if (!dir.exists(expression_out_dir)) {dir.create(expression_out_dir, recursive = TRUE)}
if (!dir.exists(sizefactor_out_dir)) {dir.create(sizefactor_out_dir, recursive = TRUE)}
if (!dir.exists(peak_info_dir)) {dir.create(peak_info_dir, recursive = TRUE)}
#------------------------------------------------------
# Count files

macs2_df <- data.frame(
  Chr = as.character(seqnames(macs2_gr)),
  Start = start(macs2_gr),
  End = end(macs2_gr),
  Geneid =macs2_gr$peakID
)

# Create ID format
macs2_df$ID <- paste(macs2_df$Chr, macs2_df$Start, macs2_df$End, sep="_")
fnf_cols <- grep("FNF", colnames(counts_macs2), value=TRUE)

# Create count table
ctl_cols <- grep("CTL", colnames(counts_macs2), value=TRUE)
ctl_count <- cbind(macs2_df, counts_macs2[, ctl_cols, drop=FALSE])

fnf_cols <- grep("FNF", colnames(counts_macs2), value=TRUE)
fnf_count <- cbind(macs2_df, counts_macs2[, fnf_cols, drop=FALSE])


# Write the full count tables for reference
write.table(ctl_count, file.path(expression_out_dir, "peaks_pbs.txt"),
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(fnf_count, file.path(expression_out_dir, "peaks_fnf.txt"),
            row.names=FALSE, quote=FALSE, sep="\t")


#------------------------------------------------------
# Create couunt file by chromosome
# Function 

save_counts_by_chromosome <- function(count_table, sample_cols, chromosomes, prefix, out_dir) {
  for (i in chromosomes) {
    message(paste0("  Processing ", i, " for ", toupper(prefix), "..."))

    # Filter by chromosome
    chr_specific <- subset(count_table, count_table$Chr == i)

    # Extract only count columns (skip the first columns with metadata)
    first_count_col <- which(colnames(chr_specific) == sample_cols[1])
    counts_data <- chr_specific[, first_count_col:ncol(chr_specific)]

    # Set row names to peak IDs
    row.names(counts_data) <- chr_specific$ID

    # Ensure all columns are numeric
    counts_data[] <- lapply(counts_data, function(x) as.numeric(x))

    # Create identifier and save matrix
    ID <- paste0(i, "_", prefix)
    countslist <- list(counts_data)
    names(countslist) <- ID

    saveRasqualMatrices(countslist, out_dir, file_suffix = "expression")
  }
}


save_counts_by_chromosome(ctl_count, ctl_cols, chromosomes, "pbs", expression_out_dir)
save_counts_by_chromosome(fnf_count, fnf_cols, chromosomes, "fnf", expression_out_dir)

#------------------------------------------------------
# Create size factor files

save_size_factors_by_chromosome <- function(norm_df, chromosomes, prefix, out_dir) {
  for (i in chromosomes) {
    message(paste0("  Processing ", i, " for ", toupper(prefix), "..."))
    
    chr_specific <- subset(norm_df, norm_df$Chr == i) # Filter by chromosome
    
    # Extract only normalization factor columns
    norm_cols <- setdiff(colnames(chr_specific), c("Chr", "PeakID")) 
    normfactors <- chr_specific[, norm_cols, drop=FALSE]

    # Set row names to peak IDs
    row.names(normfactors) <- chr_specific$PeakID

    # Ensure all columns are numeric
    normfactors[] <- lapply(normfactors, function(x) as.numeric(x))

    # Create identifier and save matrix
    ID <- paste0(i, "_", prefix)
    normlist <- list(normfactors)
    names(normlist) <- ID
    saveRasqualMatrices(normlist, out_dir, file_suffix = "size_factors")
  }
}

# Make sure CQN and the other expression cols matches
ctl_cqn_cols <- intersect(ctl_cols, colnames(cqnNormFactors))
fnf_cqn_cols <- intersect(fnf_cols, colnames(cqnNormFactors))

if (length(ctl_cqn_cols) != length(ctl_cols)) {
    warning("Not all PBS samples have normalization factors") 
} else {
  message("PBS colnames match")
}
if (length(fnf_cqn_cols) != length(fnf_cols)) {
  warning("Not all FNF samples have normalization factors")
} else {
  message("FNF colnames match")
}

CQN_pbs <- cqnNormFactors[, ctl_cqn_cols]
CQN_fnf <- cqnNormFactors[, fnf_cqn_cols]

CQN_pbs_df <- as.data.frame(CQN_pbs)
CQN_pbs_df$Chr <- macs2_df$Chr
CQN_pbs_df$PeakID <- macs2_df$ID

CQN_fnf_df <- as.data.frame(CQN_fnf)
CQN_fnf_df$Chr <- macs2_df$Chr
CQN_fnf_df$PeakID <- macs2_df$ID

save_size_factors_by_chromosome(CQN_pbs_df, chromosomes, "pbs", sizefactor_out_dir)
save_size_factors_by_chromosome(CQN_fnf_df, chromosomes, "fnf", sizefactor_out_dir)

#------------------------------------------------------
# save Peak info

save_peak_info_by_chromosome <- function(peak_df, chromosomes, out_dir, include_geneid = FALSE) {
  for (i in chromosomes) {
    message(paste0("  Creating peak info for ", i))
    chr_specific <- subset(peak_df, peak_df$Chr == i)
    if (include_geneid && "Geneid" %in% colnames(chr_specific)) {
      peak_file <- data.frame(
        GeneID = chr_specific$Geneid,
        PeakID = chr_specific$ID,
        Chr = chr_specific$Chr,
        Start = chr_specific$Start,
        End = chr_specific$End
      )
    } else {
      peak_file <- data.frame(
        PeakID = chr_specific$ID,
        Chr = chr_specific$Chr,
        Start = chr_specific$Start,
        End = chr_specific$End
      )
    }
    file_name <- file.path(out_dir, paste0("peak_info_", i, ".txt"))
    write.table(peak_file, file = file_name,
                quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}


save_peak_info_by_chromosome(macs2_df, chromosomes, peak_info_dir, include_geneid = FALSE)
#------------------------------------------------------