#This is MTC permutation correction
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(grid)  # for unit()
library(httpgd)

# Define parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
data_dir <- file.path(base_dir, "caQTL/data")
plot_dir <- file.path(base_dir, "caQTL/plots")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
#-------------------------------------------------
#Parameters
window_type <- "25kb"  # or 1,10,25, 100kb
conditions <- c("pbs", "fnf")
pc <- paste0("pc", 0)
chromosomes <- 1:22


# Function for the Alasoo's getFDR function

getFDR <- function(p1, p0, alpha=0.1, z=NULL, subset=NULL){
  if(is.null(z)){
    a=0
    for(itr in 1:10){
      a=getFDR(p1, p0, alpha, rev(a+0:100/100^itr), subset)
    }
    return(a)
  } else {
    if(!is.null(subset)){
      p1=p1[subset]
      p0=p0[subset]
    }
    p1=p1[!is.na(p1)]
    p0=p0[!is.na(p0)]
    x=NULL
    for(i in z){
      x=c(x, (sum(p0<i)/length(p0))/(sum(p1<i)/length(p1)))
    }
    return(max(c(0, z[x<alpha]), na.rm=T))
  }
}
#-------------------------------------------------
# permutation to combine them as one file

for (cond in conditions) {
    # Collect results from all chromosomes (eigenMT)
    chr_results <- lapply(chromosomes, function(chr) {
      file_path <- file.path(
        base_dir, "caQTL/eigenMT_perm/results",
        paste0("window_", window_type), cond, pc,
        sprintf("eigenMT_chr%s_%s_%s_perm1.txt", chr, cond, pc)
      )
      fread(file_path)
    })
    
    # Combine and filter out NULL results
    all_results <- rbindlist(chr_results[!sapply(chr_results, is.null)])
    all_results <- all_results[order(all_results$`p-value`),]
  
    output_dir <- file.path(
      base_dir, "caQTL/eigenMT_perm/combined_results",
      paste0("window_", window_type), cond, pc
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write combined results
    combined_file <- file.path(output_dir, sprintf("eigenMT_combined_%s_%s.txt", cond, pc))
    fwrite(all_results, combined_file, sep = "\t", quote = FALSE)
    
}

#---------------------------------------------------
# Apply getFDR
for (cond in conditions) {
  
  # Load real
  real_results <- fread(file.path(
    base_dir, "caQTL/eigenMT/combined_results",
    paste0("window_", window_type), cond, pc,
    sprintf("eigenMT_combined_%s_%s.txt", cond, pc)
  ))
  
  # Load permutation
  permutation <- fread(file.path(
    base_dir, "caQTL/eigenMT_perm/combined_results",
    paste0("window_", window_type), cond, pc,
    sprintf("eigenMT_combined_%s_%s.txt", cond, pc)
  ))
  
  # Find optimal threshold using getFDR
  # alpha = 0.05 for 5% FDR
  optimal_threshold <- getFDR(
    p1 = real_results$BF,
    p0 = permutation$BF,
    alpha = 0.05
  )
  
  # Apply to chromosomes
  for (chr in chromosomes) {
    rasqual <- fread(file.path(
      data_dir, "01_BF_adjusted",
      paste0("window_", window_type), cond, pc,
      sprintf("chr%s_%s_%s_%s_BFadjusted.txt", chr, cond, pc, window_type)
    ))
    
    rasqual_post <- rasqual[BF <= optimal_threshold]

    out_dir <- file.path(data_dir, "02_MTC_final", paste0("window_", window_type), cond)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    fwrite(rasqual_post,
           file.path(data_dir, "02_MTC_final", paste0("window_", window_type), cond,
                     sprintf("chr%s_%s_%s_MTCFinal.txt", chr, cond, pc)),
           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

# Combine
for (cond in conditions) {
  Final_DF <- rbindlist(lapply(chromosomes, function(chr) {
    fread(file.path(data_dir, "02_MTC_final", paste0("window_", window_type), cond,
                    sprintf("chr%s_%s_%s_MTCFinal.txt", chr, cond, pc)))
  }))
  
  fwrite(Final_DF,
         file.path(data_dir, "02_MTC_final", paste0("window_", window_type), cond,
                   sprintf("%s_%s_AllResults_%s_MTCFinal.txt", cond, pc, window_type)),
         sep = "\t")
}

# for (cond in conditions) {
#   permutation_path <- file.path(
#     base_dir, "caQTL/eigenMT_perm/combined_results",
#     paste0("window_", window_type), cond, pc,
#     sprintf("eigenMT_combined_%s_%s.txt", cond, pc)
#   )
  
#   permutation <- fread(permutation_path)
#   permutation <- permutation[order(permutation$`p-value`), ]
  
#   ThresholdRowCounts <- floor(bf_threshold * nrow(permutation))
#   thresholdvalue <- permutation$`p-value`[ThresholdRowCounts]
  
#   cat("  Threshold:", sprintf("%.6e", thresholdvalue), "\n")
  
#   # Apply to each chromosome
#   for (chr in chromosomes) {
    
#     # Load BF-adjusted RASQUAL for this chromosome
#     rasqual_BFpath <- file.path(
#       data_dir, "01_BF_adjusted",
#       paste0("window_", window_type), cond, pc,
#       sprintf("chr%s_%s_%s_%s_BFadjusted.txt", chr, cond, pc, window_type)
#     )
    
#     if (!file.exists(rasqual_BFpath)) {
#       cat("  Chr", chr, ": file not found\n")
#       next
#     }
    
#     rasqual <- fread(rasqual_BFpath)
    
#     # Filter by threshold
#     index <- which(rasqual$BF <= thresholdvalue)
#     rasqual_post <- rasqual[index, ]
    
#     cat("  Chr", chr, ":", nrow(rasqual_post), "significant SNPs\n")
    
#     # Save per-chromosome
#     output_dir <- file.path(
#       data_dir, "02_MTC_final",
#       paste0("window_", window_type), cond
#     )
#     dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
#     output_file <- file.path(
#       output_dir,
#       sprintf("chr%s_%s_%s_MTCFinal.txt", chr, cond, pc)
#     )
    
#     fwrite(rasqual_post, output_file, sep = "\t", 
#            row.names = FALSE, col.names = TRUE, quote = FALSE)
#   }
# }


# #---- Combine all chromosome for final output---------

# for (cond in conditions) {
#   cat("\nProcessing:", cond, "\n")
#   # Initialize combined dataframe
#   Final_DF <- NULL
  
#   for (chr in chromosomes) {    
#     filepath <- file.path(
#       data_dir, "02_MTC_final",
#       paste0("window_", window_type), cond,
#       sprintf("chr%s_%s_%s_MTCFinal.txt", chr, cond, pc)
#     )
  
#     chr_data <- fread(filepath)    
  
#     if (is.null(Final_DF)) {
#       Final_DF <- chr_data
#     } else {
#       Final_DF <- rbind(Final_DF, chr_data)
#     }
    
#     cat("  Chr", chr, ":", nrow(chr_data), "SNPs added\n")
#   }
  
#   # Save combined results
#   output_dir <- file.path(
#     data_dir, "02_MTC_final",
#     paste0("window_", window_type), cond
#   )
  
#   outputfile <- file.path(
#     output_dir,
#     sprintf("%s_%s_AllResults_25kb_MTCFinal.txt", cond, pc)
#   )
  
#   fwrite(Final_DF, outputfile, sep = "\t", 
#          row.names = FALSE, col.names = TRUE, quote = FALSE)
  
#   cat("\nTotal significant SNPs:", nrow(Final_DF), "\n")
#   cat("Unique peaks:", length(unique(Final_DF$Feature)), "\n")
#   cat("Saved to:", outputfile, "\n")
# }


# Load data
real_results <- fread(file.path(
  base_dir, "caQTL/eigenMT/combined_results",
  paste0("window_", window_type), "pbs", "pc0",
  "eigenMT_combined_pbs_pc0.txt"
))

permutation <- fread(file.path(
  base_dir, "caQTL/eigenMT_perm/combined_results",
  paste0("window_", window_type), "pbs", "pc0",
  "eigenMT_combined_pbs_pc0.txt"
))

# 1. Histogram comparison
df_plot <- rbind(
  data.frame(BF = real_results$BF, type = "Real"),
  data.frame(BF = permutation$BF, type = "Permutation")
)

p1 <- ggplot(df_plot, aes(x = BF, fill = type)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  xlim(0, 0.1) +
  scale_fill_manual(values = c("Real" = "#2057A7", "Permutation" = "gray50")) +
  labs(title = "BF Distribution: Real vs Permutation",
       x = "BF (Bonferroni corrected p-value)",
       y = "Count") +
  theme_bw()

# 2. Cumulative distribution
p2 <- ggplot(df_plot, aes(x = BF, color = type)) +
  stat_ecdf(size = 1) +
  xlim(0, 0.05) +
  scale_color_manual(values = c("Real" = "#2057A7", "Permutation" = "gray50")) +
  labs(title = "Cumulative Distribution",
       x = "BF threshold",
       y = "Proportion of peaks passing") +
  theme_bw()

# 3. FDR curve
thresholds <- seq(0, 0.05, by = 0.0001)
fdr_values <- sapply(thresholds, function(t) {
  n_real <- sum(real_results$BF <= t)
  n_perm <- sum(permutation$BF <= t)
  if (n_real == 0) return(NA)
  return(n_perm / n_real)
})

p3 <- ggplot(data.frame(threshold = thresholds, FDR = fdr_values), 
             aes(x = threshold, y = FDR)) +
  geom_line(size = 1, color = "#2057A7") +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
  labs(title = "Empirical FDR vs Threshold",
       x = "BF threshold",
       y = "Empirical FDR") +
  ylim(0, 0.3) +
  theme_bw()
