
source("/work/users/s/e/seyoun/cQTL_caQTL/scripts/atac_workflow/caQTL/response_caQTL/visQTL.r")
base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
plot_dir <- file.path(base_dir, "caQTL/plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
# Parameters
# base_dir <- "/work/users/s/e/seyoun/cQTL_caQTL/atac_output"
# response_dir <- file.path(base_dir, "caQTL/data/response_qtl")
# load(file.path(response_dir, "caQTL_prepdata.RData"))
# load(file.path(response_dir, "response_caQTL_PBS_highconf.RData"))
# load(file.path(response_dir, "response_caQTL_FNF_highconf.RData"))

# nominal rds
#pbs_rasqual <- readRDS(file.path(response_dir, "pbs_rasqual_nominal.rds"))
#fnf_rasqual <- readRDS(file.path(response_dir, "fnf_rasqual_nominal.rds"))

#------------------------------------------------------------------------------------

p_boxplot <- create_caqtl_boxplot(
  peak_id = top_pbs$peak,
  snp_id = top_pbs$snp,
  highConf = pbs_highconf
)

p_boxplot <- create_caqtl_boxplot(
  peak_id = top_fnf$peak,
  snp_id = top_fnf$snp,
  highConf = fnf_highconf
)

# ggsave("response_caQTL_boxplot_example.pdf", p_boxplot, 
#        width = 4, height = 3)
top_pbs <- pbs_highconf[5, ]


pdf(file.path(plot_dir, "highconf_pbs_all.pdf"), width = 4.5, height = 5.5)
for(i in 1:nrow(pbs_highconf)){
pageCreate(width = 4.5, height = 5.5, default.units = "inches",showGuides = FALSE)
top_pbs <- pbs_highconf[i, ]
plot_caqtl_multitrack(
  peak_id = top_pbs$peak,
  snp_id = top_pbs$snp,
  primary_dataset = "PBS",
  x_start = 0.5,
  y_start = 0.5,
  width = 3,
  height = 2.5,
  zoom_range = 50000,
  add_atac_signal = TRUE, 
  add_rna_signal = TRUE
)
}
dev.off()

pdf(file.path(plot_dir, "highconf_fnf_all.pdf"), width = 4.5, height = 5.5)
for(i in 1:nrow(fnf_highconf)){
pageCreate(width = 4.5, height = 5.5, default.units = "inches",showGuides = FALSE)
top_fnf <- fnf_highconf[i, ]
plot_caqtl_multitrack(
  peak_id = top_fnf$peak,
  snp_id = top_fnf$snp,
  primary_dataset = "FNF",
  x_start = 0.5,
  y_start = 0.5,
  width = 3.5,
  height = 2.5,
  zoom_range = 50000,
  add_atac_signal = TRUE, 
  add_rna_signal = TRUE
  )
}
dev.off()


pdf(file.path(plot_dir, "poster_fig_case.pdf"), width = 6.5, height = 5.5)
pageCreate(width = 6.5, height = 5.5, default.units = "inches",showGuides = FALSE)
top_pbs <- pbs_highconf[1, ]
plot_caqtl_multitrack(
  peak_id = top_pbs$peak,
  snp_id = top_pbs$snp,
  primary_dataset = "PBS",
  x_start = 0.5,
  y_start = 0.5,
  width = 2.5,
  height = 1.4,
  zoom_range = 50000,
  add_atac_signal = TRUE, 
  add_rna_signal = TRUE
)

#I need to change the line # 227-228 to:
#minregion <- max(1, peak_center - zoom_range  ) # max(1, peak_center - zoom_range  - zoom_range/2
#maxregion <- peak_center + zoom_range # peak_center + zoom_range/2

top_fnf <- fnf_highconf[11, ]
plot_caqtl_multitrack(
  peak_id = top_fnf$peak,
  snp_id = top_fnf$snp,
  primary_dataset = "FNF",
  x_start = 3.7,
  y_start = 0.5,
  width = 2.5,
  height = 1.4,
  zoom_range = 50000,
  add_atac_signal = TRUE, 
  add_rna_signal = TRUE
)

# pbs_boxplot <- create_caqtl_boxplot(
#   peak_id = top_pbs$peak,
#   snp_id = top_pbs$snp,
#   highConf = pbs_highconf
# )
# save(pbs_boxplot, file = file.path(plot_dir, "pbs_caqtl_case_boxplot.rds"))
# fnf_boxplot <- create_caqtl_boxplot(
#   peak_id = top_fnf$peak,
#   snp_id = top_fnf$snp,
#   highConf = fnf_highconf
# )
# save(fnf_boxplot, file = file.path(plot_dir, "fnf_caqtl_case_boxplot.rds"))

plotGG(pbs_boxplot, x = 0.1, y = 4, width = 3, height = 1.25)
plotGG(fnf_boxplot, x = 3.3, y = 4, width = 3, height = 1.25)
dev.off()
