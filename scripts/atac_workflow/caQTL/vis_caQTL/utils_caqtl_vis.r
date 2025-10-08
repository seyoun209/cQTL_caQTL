# Function

# Load libraries ----------------------------------------------------------
# Find the function to caQTL LD variants and combine. 

get_caqtl_ld_variants <- function(caPeaks, ld_dir, r2_threshold = 0.5) {
  # caPeaks: data.table with caQTL info (must include column 'snp')
  # Returns: data.table with caQTL SNP, LD variants, and R2 values
  
  result_list <- list()
  
  for(i in seq_len(nrow(caPeaks))) {
    snp_id <- caPeaks$snp[i]
    
    # Extract chromosome
    chr_sub <- strsplit(snp_id, ":")[[1]][1]
    ld_file <- file.path(ld_dir, chr_sub, paste0(snp_id, ".ld"))
    
    if(file.exists(ld_file)) {
      tryCatch({
        ld_data <- fread(ld_file)
        
        # Filter by R2 threshold (0.5 instead of 0.7)
        ld_hi <- ld_data[R2 > r2_threshold]
        
        if(nrow(ld_hi) > 0) {
          # Add the lead SNP info
          ld_hi$lead_snp <- snp_id
          ld_hi$peak <- caPeaks$peak[i]
          result_list[[i]] <- ld_hi
        }
      }, error = function(e) {
        message(paste("Error processing", snp_id, ":", e$message))
      })
    }
  }
  
  # Combine all results
  if(length(result_list) > 0) {
    return(rbindlist(result_list, fill = TRUE))
  } else {
    return(data.table())
  }
}

