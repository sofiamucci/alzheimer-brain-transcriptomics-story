library(DESeq2)
library(plotly)
library(dplyr)

## Load data from 05_qc.R
vsd_data <- readRDS("data/vsd_data.rds")
vsd_mat_data <- assay(vsd_data)
col.data <- readRDS("data/col_data.rds")

## Full PCA across all components
pca_full <- prcomp(t(vsd_mat_data), scale = TRUE)
summary(pca_full)$importance[, 1:6]

# Note: PC1-PC2 alone do not show clean separation between conditions 
# (only ~36% of variance explained).

## Interactive 3D PCA (PC1, PC2, PC3)
pca_df <- data.frame(pca_full$x[, 1:5], condition = col.data$condition)

plot_ly(pca_df, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~condition, 
        colors = c("Young" = "darkturquoise", "Old" = "grey40", "AD" = "darkorange"),
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 6)) %>%
  layout(title = "PCA: PC1, PC2, PC3 (interactive)")


plot_ly(pca_df, 
        type = "splom",
        dimensions = list(
          list(label = "PC1", values = ~PC1),
          list(label = "PC2", values = ~PC2),
          list(label = "PC3", values = ~PC3),
          list(label = "PC4", values = ~PC4),
          list(label = "PC5", values = ~PC5)
        ),
        color = ~condition,
        colors = c("Young" = "darkturquoise", "Old" = "grey40", "AD" = "darkorange"),
        marker = list(size = 5)) %>%
  layout(title = "PCA: PC1 to PC5 (interactive matrix)")

## Static pairs plot (PC1-PC5)
# Color palette for the 3 conditions
colors <- c("Young" = "darkturquoise", "Old" = "grey40", "AD" = "darkorange")

par(oma = c(0, 0, 3, 0))

pairs(pca_df[, 1:5], 
      col = colors[pca_df$condition], 
      pch = 19, 
      cex = 1.2,
      main = "PCA: PC1 to PC5")

mtext("turquoise = Young, grey = Old, orange = AD", 
      side = 3, line = 0.1, cex = 0.8, outer = TRUE)

# Reset margins for future plots
par(oma = c(0, 0, 0, 0))











  

  