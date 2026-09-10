library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(pROC)

## Load metadata and salmon/tximport quantifications
col.data <- read.csv("data/sample_metadata_full.csv")
txi <- readRDS("data/txi_full.rds")

# Keep the original SRR accession as a column before renaming
col.data$srr_accession <- col.data$external_id

# Build readable sample names: condition + replicate number within each group
col.data <- col.data %>%
  group_by(condition) %>%
  mutate(sample_name = paste0(condition, "_", row_number())) %>%
  ungroup() %>%
  as.data.frame()

rownames(col.data) <- col.data$sample_name
colnames(txi$counts) <- col.data$sample_name[match(colnames(txi$counts), col.data$srr_accession)]
colnames(txi$abundance) <- colnames(txi$counts)
colnames(txi$length) <- colnames(txi$counts)

# sanity check: samples match and are in the same order
all(colnames(txi$counts) == rownames(col.data))
col.data[, c("srr_accession", "condition", "sample_name")]

# counts_data: raw counts matrix extracted from tximport
counts_data <- txi$counts
counts_data <- round(counts_data)

## Build the DESeq2 dataset object directly from tximport
data.dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                   colData = col.data,
                                   design = ~ condition)
data.dds

# Filter out low-count genes
keep <- rowSums(counts(data.dds)) >= 10
data.ddsk <- data.dds[keep, ]

# Set reference level (Young = baseline for aging/disease comparisons)
data.ddsk$condition <- relevel(data.ddsk$condition, ref = "Young")

# Run DESeq2
data.analisis <- DESeq(data.ddsk)
resultsNames(data.analisis)

# Old vs AD: the disease-specific comparison, controlling for baseline aging
res_old_vs_ad <- results(data.analisis, contrast = c("condition", "AD", "Old"))
res_old_vs_ad

plotDispEsts(data.analisis)
summary(res_old_vs_ad)

# Normalize and estimate size factors
estimated_data.dds <- estimateSizeFactors(data.ddsk)
sizeFactors(estimated_data.dds)
normalized_data <- counts(estimated_data.dds, normalized = TRUE)

# Save normalized counts for downstream scripts
saveRDS(normalized_data, "data/normalized_data.rds")
write.csv(normalized_data, "data/normalized_data.csv", row.names = TRUE)

# Variance stabilizing transformation
vsd_data <- vst(estimated_data.dds, blind = TRUE)
vsd_mat_data <- assay(vsd_data)

# Save vsd object for use in 06_pca_multidim_exploration.R
saveRDS(vsd_data, "data/vsd_data.rds")

# Sample-to-sample correlation
vsd_cor_data <- cor(vsd_mat_data)
View(vsd_cor_data)

### PCA plot 2D
pca_data <- plotPCA(vsd_data, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

# Scree plot
PCA_2 <- prcomp(t(vsd_mat_data), scale = TRUE)
names(PCA_2)

prop_variance <- PCA_2$sdev^2 / sum(PCA_2$sdev^2)
prop_variance

# Screeplot max=1
ggplot(data = data.frame(prop_variance, pc = 1:length(prop_variance)), aes(x = pc, y = prop_variance)) +
  geom_col(fill = "red", color = "red", width = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Principal components", y = "Proportion of variance explained")

# Screeplot max=0.25
ggplot(data = data.frame(prop_variance, pc = 1:length(prop_variance)), aes(x = pc, y = prop_variance)) +
  geom_col(fill = "red", color = "red", width = 0.5) +
  scale_y_continuous(limits = c(0, 0.25)) +
  labs(x = "Principal components", y = "Proportion of variance explained")

# Sample-to-sample correlation heatmap
pheatmap(vsd_cor_data)

# Save DESeq2 results object for use in 07_volcano_plot.R
saveRDS(res_old_vs_ad, "data/res_data.rds")

# Save data for use in downstream scripts
saveRDS(col.data, "data/col_data.rds")
