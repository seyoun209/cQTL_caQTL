# Combine all chromosome files into one file for each condition
library(data.table)
library(dplyr)

base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
data_dir <- file.path(base_dir, "caQTL/data")
window_type <- "25kb"
conds <- c("pbs", "fnf")
chromosomes <- 1:22

for (cond in conds) {  
  # Directory containing chromosome files
  cond_dir <- file.path(data_dir, "05_LDFiltered_LeadsOnly", paste0("window_", window_type))
  
  # Collect all existing files
  files <- sprintf("%s/%s_chr%s_LDfiltered_leads_R2lt02.csv", cond_dir, cond, chromosomes)
  files <- files[file.exists(files)]
    
  # Read all chromosome files
  merged <- rbindlist(lapply(files, fread), fill = TRUE)
  
  # Drop placeholder rows if any
  merged <- merged %>% filter(!is.na(Feature))
  
  # Save combined output
  out_file <- file.path(cond_dir, sprintf("%s_AllChr_LDfiltered_leads_R2lt02.csv", cond))
  fwrite(merged, out_file)
}

#-------------------------------------------------------------------
# Save the lead snps -----------------------------------------------

#PBS
cond_dir <- file.path(data_dir, "05_LDFiltered_LeadsOnly", paste0("window_", window_type))
  
# Collect all existing files
pbs_snps <- fread(file.path(cond_dir, sprintf("%s_AllChr_LDfiltered_leads_R2lt02.csv", "pbs")))
pbs_leadSnps <- pbs_snps %>%
  group_by(Feature) %>%
  arrange(BF, Distance_from_PeakCenter) %>%
  slice(1) %>%
  ungroup()

fnf_snps <- fread(file.path(cond_dir, sprintf("%s_AllChr_LDfiltered_leads_R2lt02.csv", "fnf")))
fnf_leadSnps <- fnf_snps %>%
  group_by(Feature) %>%
  arrange(BF, Distance_from_PeakCenter) %>%
  slice(1) %>%
  ungroup()



