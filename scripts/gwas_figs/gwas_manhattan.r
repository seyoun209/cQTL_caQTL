library(plotgardener)
library(stringr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggtext)
library(colorspace)
library(RColorBrewer)


# ---------- parameters ----------
ss_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/external/Hatzikotoulas_2025/hg38"
oa_subtypes <- c("ALLOA", "FINGER", "HAND", "HIP", "HIPKNEE", "KNEE", "SPINE", "THR", "THUMB", "TJR", "TKR")


read_gwas_subtype <- function(subtype, dir = ss_dir){
  subdir <- file.path(dir, subtype, "summary_statistics")
  search_dir <- if (dir.exists(subdir)) subdir else dir
  
  files <- list.files(search_dir, pattern = paste0("^", subtype, "_chr.*liftOver\\.csv$"), full.names = TRUE)
  if(length(files) == 0) return(NULL)
  dt_list <- lapply(files, function(f){
    dt <- fread(f)
    dt[, pos := as.integer(hg38pos)]
    dt[, chrom := ifelse(grepl("^chr", CHR, ignore.case = TRUE), CHR, paste0("chr", CHR))]
    dt[, p := as.numeric(P)]
    dt[p == 0, p := 1e-15]
    dt[, snp := as.character(CPTID)]
    dt <- dt[!is.na(p) & !is.na(pos) & !is.na(chrom)] |> select(chrom, pos, p, snp)
    dt
  })
  dt_list <- dt_list[!sapply(dt_list, is.null)]
  if(length(dt_list) == 0) return(NULL)
  dt_all <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  dt_all
}



# ---------- single-subtype manhattan example (ALLOA) ----------
# dt_alloa <- read_gwas_subtype("ALLOA")

# gwas_allOA <- data.frame(chrom = dt_alloa$chrom,
#                     pos   = dt_alloa$pos,
#                     p     = dt_alloa$p,
#                     snp   = dt_alloa$snp,
#                     stringsAsFactors = FALSE)
# # simple LD placeholder
# gwas_allOA$LD <-  NA_real_

# if (nrow(gwas_allOA) > 200000) {
#   top  <- gwas_allOA %>% arrange(p) %>% slice_head(n = 5000)
#   samp <- gwas_allOA %>% sample_n(95000)
#   gwas_plot <- bind_rows(top, samp)
# } else {
#   gwas_plot <- gwas_allOA
# }

# pageCreate(width = 7.5, height = 4.5, default.units = "inches", showGuides = FALSE)

# manhattan_all <- plotManhattan(
#     data = gwas_plot,
#     assembly = "hg38",
#     fill = c("grey", "#37a7db"),
#     sigLine = TRUE,
#     trans = "-log10",
#     col = "grey", lty = 2, range = c(0, 14),
#     x = 0.5, y = 0, width = 6.5, height = 2,
#     just = c("left", "top"),
#     default.units = "inches"
# )
# annoGenomeLabel(
#     plot = manhattan_all, x = 0.5, y = 2, fontsize = 8,
#     just = c("left", "top"),
#     default.units = "inches"
# )

# plotText(
#     label = "Chromosome", fontsize = 8,
#     x = 3.75, y = 2.20, just = "center", default.units = "inches"
# )

# ## Annotate y-axis
#     annoYaxis(
#         plot = manhattan_all, at = c(0, 2, 4, 6, 8, 10, 12, 14,16,18,20,22,24,26,28,30),
#         axisLine = TRUE, fontsize = 8
#     )

# ## Plot y-axis label
# plotText(
#     label = "-log10(p-value)", x = 0.15, y = 1, rot = 90,
#     fontsize = 8, fontface = "bold", just = "center",
#     default.units = "inches"
# )


out_dir <- file.path("/work/users/s/e/seyoun/cQTL_caQTL/atac_output/plots")

for(sub in oa_subtypes){
    print(sub)
    dt_gwas <- read_gwas_subtype(sub)
    if (is.null(dt_gwas) || nrow(dt_gwas) == 0) next
    gwas <- data.frame(chrom = dt_gwas$chrom, pos = dt_gwas$pos, p = dt_gwas$p, snp = dt_gwas$snp, stringsAsFactors = FALSE)
    gwas$LD <-  NA_real_
out_png <- file.path(out_dir, paste0(sub, "_gwas_hg38_summarystat.png"))
png(filename = out_png, width = 7.5 * 300, height =  4.5 * 300, res = 300,  bg = "transparent")
pageCreate(width = 7.5, height = 4.5, default.units = "inches", showGuides = FALSE)
manhattan_all <- plotManhattan(
    data = gwas,
    assembly = "hg38",
    fill = c("grey", "#37a7db"),
    sigLine = TRUE,
    trans = "-log10",
    col = "grey", lty = 2, range = c(0, 30),
    x = 0.5, y = 0.5, width = 6.5, height = 2,
    just = c("left", "top"),
    default.units = "inches"
)


annoGenomeLabel(
    plot = manhattan_all, x = 0.5, y = 2.5, fontsize = 8,
    just = c("left", "top"),
    default.units = "inches"
)

plotText(
    label = "Chromosome", fontsize = 8,
    x = 3.75, y = 2.70, just = "center", default.units = "inches"
)



## Annotate y-axis
    annoYaxis(
        plot = manhattan_all, at = c(0, 2, 4, 6, 8, 10, 12, 14,16,18,20,22,24,26,28,30),
        axisLine = TRUE, fontsize = 8
    )

## Plot y-axis label
plotText(
    label = "-log10(p-value)", x = 0.15, y = 1.5, rot = 90,
    fontsize = 8, fontface = "bold", just = "center",
    default.units = "inches"
)  
dev.off()
}

