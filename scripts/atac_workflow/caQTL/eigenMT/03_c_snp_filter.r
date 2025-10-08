library(data.table)
library(dplyr)
library(stringr)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 03_c_LDsnpFilter.r <chromosome> <condition>")
chr <- as.integer(args[1])
cond <- args[2]


# Define parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
data_dir <- file.path(base_dir, "caQTL/data")
plot_dir <- file.path(base_dir, "caQTL/plots")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------
# Parameters
window_type <- "25kb"
#conditions <- c("pbs", "fnf")
#conditions <- c("fnf")
pc <- "pc0"
#chromosomes <- 1:22

out_dir <- file.path(data_dir, "05_LDFiltered_LeadsOnly", paste0("window_", window_type))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


rasqual_file <- file.path(
  data_dir, "02_MTC_final", paste0("window_", window_type), cond,
  sprintf("%s_%s_AllResults_%s_MTCFinal.txt", cond, pc, window_type)
)
cat("Loading:", rasqual_file, "\n")

RASQUAL <- fread(rasqual_file)
RASQUAL <- RASQUAL %>%
  mutate(Distance_from_PeakCenter = abs(Distance_From_Peak))

CHR <- paste0("chr", chr)
RASQUAL_chr <- RASQUAL %>% filter(Chromosome == CHR)
if (nrow(RASQUAL_chr) == 0) {
  cat("No data for", CHR, "\n")
  quit(save = "no")
}

#-------------------------------------------------
# Load LD file
#-------------------------------------------------
ld_file <- file.path(
  data_dir, "04_LD_results", paste0("window_", window_type), "ld0",
  sprintf("chr%s.LD.ld", chr)
)
cat("Loading LD:", ld_file, "\n")

LD <- fread(ld_file)
names(LD) <- make.names(names(LD), unique = TRUE)

# If first row repeats 
if (all(names(LD) == as.character(unlist(LD[1, ])))) {
  LD <- LD[-1, ]
}


# Pre-filter LD
snps_for_filter <- RASQUAL_chr$rs_ID
LD_main <- LD %>% filter(SNP_A %in% snps_for_filter & SNP_B %in% snps_for_filter)

#-------------------------------------------------
# LD filtering per peak
#-------------------------------------------------
Final_LDFiltered <- data.frame()
peaks <- unique(RASQUAL_chr$Feature)

for (peak in peaks) {
  print(peak)
  peak_data <- RASQUAL_chr %>% filter(Feature == peak)
  if (nrow(peak_data) == 0) next
  peak_data <- peak_data[order(peak_data$BF, peak_data$Distance_from_PeakCenter), ]
  Final_LDFiltered <- rbind(Final_LDFiltered, peak_data[1, ])

  if (nrow(peak_data) > 1) {
    lead <- peak_data$rs_ID[1]
    SecondarySNPs <- peak_data$rs_ID[-1]

    LD_lead <- LD_main %>%
      filter((SNP_A == lead & SNP_B %in% SecondarySNPs) |
             (SNP_B == lead & SNP_A %in% SecondarySNPs)) %>%
      filter(R2 < 0.2)

    while (nrow(LD_lead) > 0) {
      snp_list <- unique(c(LD_lead$SNP_A, LD_lead$SNP_B))
      sub <- peak_data %>% filter(rs_ID %in% snp_list & rs_ID != lead)
      if (nrow(sub) == 0) break
      sub <- sub[order(sub$BF, sub$Distance_from_PeakCenter), ]
      Final_LDFiltered <- rbind(Final_LDFiltered, sub[1, ])

      lead <- sub$rs_ID[1]
      SecondarySNPs <- sub$rs_ID[-1]
      LD_lead <- LD_main %>%
        filter((SNP_A == lead & SNP_B %in% SecondarySNPs) |
               (SNP_B == lead & SNP_A %in% SecondarySNPs)) %>%
        filter(R2 < 0.2)
    }
  }
}

Final_LDFiltered <- Final_LDFiltered %>% filter(!is.na(Feature))
out_file <- file.path(out_dir, sprintf("%s_chr%s_LDfiltered_leads_R2lt02.csv", cond, chr))
fwrite(Final_LDFiltered, out_file)
cat("Saved:", out_file, "\n")
