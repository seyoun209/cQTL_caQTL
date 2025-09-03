#!/usr/bin/env Rscript
## ATACseq QC Script
## Updated: Sept 2025
##------------------------------------------------------------------------------

# Load required libraries with error handling
required_packages <- c("ATACseqQC", "BSgenome.Hsapiens.UCSC.hg38", 
                      "TxDb.Hsapiens.UCSC.hg38.knownGene", "ChIPpeakAnno", 
                      "Rsamtools")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Required package", pkg, "is not installed"))
  }
}

# Parse command line arguments
args <- commandArgs(TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript ATACseq_QC.R <bam_file> <output_directory> <name_suffix>")
}

bam_file <- args[1]      # Input BAM file path (full path)
output_dir <- args[2]    # Output directory (full path)
name_suffix <- args[3]   # Suffix to remove from filename (e.g., ".sorted_final.bam")

cat("Input BAM file:", bam_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Name suffix:", name_suffix, "\n")

# Verify input file exists
if (!file.exists(bam_file)) {
  stop(paste("Input BAM file does not exist:", bam_file))
}

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Create output label from filename by removing the suffix
base_name <- basename(bam_file)
bam_label <- gsub(name_suffix, "", base_name)

cat("Processing BAM file:", bam_file, "\n")
cat("Output label will be:", bam_label, "\n")

# Run BAM QC analysis
cat("Running ATACseqQC bamQC analysis...\n")

tryCatch({
  # Run bamQC function
  bamqc_result <- bamQC(bam_file, outPath = NULL)
  
  # Save results
  output_file <- file.path(output_dir, paste0(bam_label, "_bamQC"))
  save(bamqc_result, file = output_file)
  
  cat("bamQC analysis completed successfully!\n")
  cat("Results saved to:", output_file, "\n")
  
  # Print summary statistics if available
  if (exists("bamqc_result") && !is.null(bamqc_result)) {
    cat("BAM QC Summary:\n")
    print(summary(bamqc_result))
  }
  
}, error = function(e) {
  cat("Error during bamQC analysis:", e$message, "\n")
  quit(status = 1)
})

cat("ATACseq QC analysis completed for sample:", bam_label, "\n")