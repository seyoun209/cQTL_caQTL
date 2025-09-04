library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(e1071)

setwd("/work/users/s/e/seyoun/CQTL_AI/Genopipe/output_biallelic/ancestry/")
## Read in pca data
pcaData <- fread("./CQTL_COA9_10_ref_merge.pca.evec", data.table = FALSE)
## Grab sample names and first 2 PC's
pcaData <- pcaData[,c(1:6)]
colnames(pcaData) <- c("sample", "PC1","PC2","PC3","PC4","PC5")

## Read in panel data
panel <- fread("/proj/phanstiel_lab/References/genomes/1000G/GRCh38/igsr-1000 genomes on grch38.tsv", data.table = FALSE) |>
  select(sample ="Sample name",
         pop ="Population code",
         super_pop ="Superpopulation code",
         gender ="Sex") 
## Match ref panel to pca data based on "sample" column

pca_panel <- left_join(pcaData, panel, by = c("sample"))


## SVN to infer super_pop
pca_panel_train <- pca_panel %>% filter(!is.na(super_pop))
pca_panel_train$super_pop <- as.factor(pca_panel_train$super_pop)
pca_panel_test <- pca_panel %>% filter(is.na(super_pop))
svm_ancestry <- svm(super_pop~PC1+PC2+PC3+PC4, data = pca_panel_train,
                    type = "C-classification", kernel = "radial")



prediction <- predict(svm_ancestry, pca_panel_test[,c("PC1", "PC2","PC3","PC4")])
pca_panel_test$super_pop <- prediction

pca_panel_test %>% dplyr::select(sample, super_pop) %>%
  dplyr::rename("Donor" = sample) %>%
  dplyr::rename("Predicted_Ancestry" = super_pop) %>%
  write_csv(file = "Current_study_predictedAncestry.csv")


## Rename our population name to argument name
pca_panel[which(is.na(pca_panel$super_pop)), "super_pop"] <- "Current_study"

## Separate panel data from sample data to control order of plotting points
panel_df <- pca_panel[which(pca_panel$super_pop != "Current_study"),] |> filter(!super_pop == "EUR,AFR")
pca_df <- pca_panel[which(pca_panel$super_pop == "Current_study"),]

## Change super_pop columns to factors to color by values
panel_df$super_pop <- as.factor(panel_df$super_pop)
pca_df$super_pop <- as.factor(pca_df$super_pop)

## Plot PC1 vs PC2, coloring by super population
pcaplot <- ggplot(panel_df, aes(PC2, PC3)) + 
  geom_point(aes(color = factor(super_pop))) + 
  geom_point(data = pca_df, aes(PC2, PC3, fill = factor(super_pop))) +
  theme_light() + labs(color = "Population", fill = NULL)

# ggplot(panel_df, aes(PC1, PC2)) + 
#   geom_point(aes(color = pop)) +
#   scale_color_manual(values = c(
#     "ACB" = "black", "ASW" = "black", "BEB" = "black", 
#     "CDX" = "black", "CEU" = "black", "CHB" = "black",
#     "CHS" = "black", "CLM" = "black", "ESN" = "black", 
#     "FIN" = "black", "GBR" = "black", "GIH" = "black",
#     "GWD" = "black", "IBS" = "black", "ITU" = "black",
#     "JPT" = "black", "KHV" = "black", "LWK" = "black",
#     "MSL" = "black", "MXL" = "black", "PEL" = "black",
#     "PJL" = "black", "PUR" = "grey", "STU" = "black",
#     "TSI" = "black", "YRI" = "black"
#   )) +
#   geom_point(data = pca_df, aes(PC1, PC2, fill = factor(super_pop))) +
#   theme_light() + labs(color = "Population", fill = NULL)
ggsave(filename = "Current_study_ancestry.pdf", plot = pcaplot)


eval_data <- fread("CQTL_COA9_10_ref_merge.eval",data.table = FALSE)
# Get first 6 eigenvalues
eval_data_subset <- head(eval_data$V1, 10)

# Create data frame for plotting
eval_df <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  Value = eval_data_subset
)

# Create scree plot
ggplot(eval_df, aes(x = PC, y = Value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(group = 1)) +
  geom_point() +
  labs(x = "Principal Components",
       y = "Eigenvalue",
       title = "Scree Plot") +
  theme_minimal()

# Using UMAP for visualization
library(umap)

# Combine PCs 1-5 from both datasets
combined_pcs <- rbind(
  panel_df[, c("PC1", "PC2", "PC3","PC4")],
  pca_df[, c("PC1", "PC2", "PC3","PC4")]
)

# Run UMAP
umap_result <- umap(combined_pcs)

# Add UMAP coordinates back to original data
panel_df$UMAP1 <- umap_result$layout[1:nrow(panel_df), 1]
panel_df$UMAP2 <- umap_result$layout[1:nrow(panel_df), 2]
pca_df$UMAP1 <- umap_result$layout[(nrow(panel_df)+1):nrow(umap_result$layout), 1]
pca_df$UMAP2 <- umap_result$layout[(nrow(panel_df)+1):nrow(umap_result$layout), 2]

# Plot UMAP results
ggplot(panel_df, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = factor(super_pop))) +
  geom_point(data = pca_df, aes(fill = factor(super_pop))) +
  theme_light() +
  labs(color = "Population", fill = NULL)


library(plotly)
plot_ly() %>%
  add_trace(
    data = panel_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers",
    color = ~super_pop,  # This already uses super_pop values
    name = ~super_pop    # Changed from "Reference Panel" to use super_pop values
  ) %>%
  add_trace(
    data = pca_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers",
    marker = list(size = 5),
    color = ~super_pop,
    name = ~super_pop    # Changed from "Current Study" to use super_pop values
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ),
    showlegend = TRUE
  )



# Create a function for multiple PC plots
plot_pc_pairs <- function(panel_df, pca_df, pc_pairs) {
  plots <- list()
  for(pair in pc_pairs) {
    p <- ggplot(panel_df, aes(.data[[pair[1]]], .data[[pair[2]]])) + 
      geom_point(aes(color = factor(super_pop))) +
      geom_point(data = pca_df, aes(fill = factor(super_pop))) +
      theme_light() +
      labs(color = "Population", fill = NULL,
           x = pair[1], y = pair[2])
    plots[[paste(pair[1], pair[2])]] <- p
  }
  return(plots)
}

# Create plots for different PC combinations
pc_pairs <- list(
  c("PC1", "PC2"),
  c("PC2", "PC3"),
  c("PC3", "PC4")
)

# Generate plots
plots <- plot_pc_pairs(panel_df, pca_df, pc_pairs)

# Arrange plots in a grid
library(gridExtra)
grid.arrange(grobs = plots, ncol = 2)
