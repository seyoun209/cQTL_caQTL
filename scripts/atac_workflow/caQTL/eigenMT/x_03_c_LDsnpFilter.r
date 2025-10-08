library(data.table)
library(dplyr)
library(stringr)
library(vroom)


# Define parameters
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
data_dir <- file.path(base_dir, "caQTL/data")
plot_dir <- file.path(base_dir, "caQTL/plots")
out_dir <- file.path(data_dir, "05_LDFiltered_LeadsOnly", paste0("window_", window_type))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------
# Parameters
window_type <- "25kb"
conditions <- c("pbs", "fnf")
#conditions <- c("fnf")
pc <- "pc0"
#chromosomes <- 1:22


#-------------------------------------------------
# Load LD snp and filtering

for (cond in conditions) {
  cat("\nProcessing:", cond, "\n")
  # Load MTC results
  rasqual_file <- file.path(
    data_dir, "02_MTC_final", paste0("window_", window_type), cond,
    sprintf("%s_%s_AllResults_25kb_MTCFinal.txt", cond, pc)
  )
  
  RASQUAL <- fread(rasqual_file)
  
  # Add distance from peak center
  RASQUAL <- RASQUAL %>%
    mutate(Distance_from_PeakCenter = abs(Distance_From_Peak))
  
  Final_LDFiltered <- data.frame()
  
  for (chr in chromosomes) {
    cat("  Chromosome:", chr, "\n")

    CHR <- paste0("chr", chr)
    RASQUAL_chr <- RASQUAL %>% filter(Chromosome == CHR)

    # Load LD file
    ld_file <- file.path(
      data_dir, "04_LD_results", paste0("window_", window_type), "ld0",
      sprintf("chr%s.LD.ld", chr)
    )

    LD <- fread(ld_file)

    # Pre-filter LD to only relevant SNPs
    snps_for_filter <- RASQUAL_chr$rs_ID
    LD_main <- LD %>% filter(SNP_A %in% snps_for_filter & SNP_B %in% snps_for_filter)

    # Loop over peaks
    peaks <- unique(RASQUAL_chr$Feature)

    for (peak in peaks) {
      peak_data <- RASQUAL_chr %>% filter(Feature == peak)
      peak_data <- peak_data[order(peak_data$BF, peak_data$Distance_From_Peak), ]
      Final_LDFiltered <- rbind(Final_LDFiltered, peak_data[1, ])

      if (nrow(peak_data) > 1) {
        lead <- peak_data$rs_ID[1]
        SecondarySNPs <- peak_data$rs_ID[-1]

        LD_lead <- LD_main %>%
          filter((SNP_A == lead & SNP_B %in% SecondarySNPs) | (SNP_B == lead & SNP_A %in% SecondarySNPs))

        LD_lead <- LD_lead %>% filter(R2 < 0.2)

        while (nrow(LD_lead) > 0) {
          snp_list <- unique(c(LD_lead$SNP_A, LD_lead$SNP_B))
          sub <- peak_data %>% filter(rs_ID %in% snp_list & rs_ID != lead)
          if (nrow(sub) == 0) break

          sub <- sub[order(sub$BF, sub$Distance_From_Peak), ]
          Final_LDFiltered <- rbind(Final_LDFiltered, sub[1, ])

          lead <- sub$rs_ID[1]
          SecondarySNPs <- sub$rs_ID[-1]

          LD_lead <- LD_main %>%
            filter((SNP_A == lead & SNP_B %in% SecondarySNPs) | (SNP_B == lead & SNP_A %in% SecondarySNPs)) %>%
            filter(R2 < 0.2)
        }
      }
    } 
  }

  # Remove placeholders or NA rows
  Final_LDFiltered <- Final_LDFiltered %>% filter(!is.na(Feature))

  # Save
  out_file <- file.path(out_dir, sprintf("%s_LDfiltered_leads_R2lt02.csv", cond))
  fwrite(Final_LDFiltered, out_file)
  cat("Saved:", out_file, "\n")
}
