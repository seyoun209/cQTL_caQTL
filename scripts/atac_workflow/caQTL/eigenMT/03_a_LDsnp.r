# Make LD SNPs list and find the LD 
library(data.table)
library(tidyverse)


# Define parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
data_dir <- file.path(base_dir, "caQTL/data")
plot_dir <- file.path(base_dir, "caQTL/plots")


dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
#-------------------------------------------------
#Parameters
window_type <- "25kb"  # or 1, 10, 100, 25kb
conditions <- c("pbs", "fnf")
pc <- paste0("pc", 0)
chromosomes <- 1:22


snp_dir <- file.path(data_dir, "03_LD_snplist", paste0("window_", window_type))
dir.create(snp_dir, recursive = TRUE, showWarnings = FALSE)
#-------------------------------------------------

#------ Make all the LD list of snps---------------
Full_RASQUAL <- NULL
for (cond in conditions) {
  file <- file.path(data_dir, "02_MTC_final", 
                    paste0("window_", window_type), cond,
                    sprintf("%s_%s_AllResults_%s_MTCFinal.txt", cond, pc, window_type))
  data <- fread(file)
  Full_RASQUAL <- rbind(Full_RASQUAL, data)
}

#------- Make SNP list per chromosome and save them----------------

for (chr in chromosomes) {
  chr_data <- Full_RASQUAL[Chromosome == paste0("chr", chr)]
  snp_list <- unique(chr_data$rs_ID)

  output_file <- file.path(snp_dir, sprintf("chr%s_SNPList.txt", chr))
  write.table(snp_list, output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Chr", chr, ":", length(snp_list), "unique SNPs\n")
}

