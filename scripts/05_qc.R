library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(pROC)
library(apeglm)

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

## Boxplot of log2(counts) distribution per sample
# Confirms that overall transcript abundance distributions are 
# comparable across samples before proceeding with normalization/analysis
logcounts <- log2(counts(data.ddsk) + 1)

condition_colors <- c("Young" = "darkturquoise", "Old" = "grey40", "AD" = "darkorange")
sample_colors <- condition_colors[col.data$condition]

par(oma = c(0, 0, 4, 0))

boxplot(logcounts, xlab = "Samples", ylab = "Log2(Counts + 1)", 
        las = 2, col = sample_colors, cex.axis = 0.6)

mtext("Relative abundance of counts", side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)

legend(x = "top", inset = -0.15, legend = names(condition_colors), 
       fill = condition_colors, horiz = TRUE, xpd = NA, bty = "n", cex = 0.9)

par(oma = c(0, 0, 0, 0))


# Set reference level (Young = baseline for aging/disease comparisons)
data.ddsk$condition <- relevel(data.ddsk$condition, ref = "Young")

# Run DESeq2
data.analisis <- DESeq(data.ddsk)
resultsNames(data.analisis)

# Old vs AD: the disease-specific comparison, controlling for baseline aging
res_old_vs_ad <- results(data.analisis, 
                         contrast = c("condition", "AD", "Old"), 
                         alpha = 0.05)
metadata(res_old_vs_ad)$alpha
summary(res_old_vs_ad)

# Old vs Young: pure aging effect (no disease)
res_old_vs_young <- results(data.analisis, 
                            contrast = c("condition", "Old", "Young"), 
                            alpha = 0.05)
metadata(res_old_vs_young)$alpha
summary(res_old_vs_young)

# AD vs Young: combined aging + disease effect
res_ad_vs_young <- results(data.analisis, 
                           contrast = c("condition", "AD", "Young"), 
                           alpha = 0.05)
metadata(res_ad_vs_young)$alpha
summary(res_ad_vs_young)

plotDispEsts(data.analisis)

# Shrink log2FoldChange estimates (type="normal" for consistency 
# across all three contrasts, avoiding apeglm/ashr dependency issues)
res_old_vs_young <- lfcShrink(data.analisis, coef = "condition_Old_vs_Young", 
                              res = res_old_vs_young, type = "normal")
res_ad_vs_young <- lfcShrink(data.analisis, coef = "condition_AD_vs_Young", 
                             res = res_ad_vs_young, type = "normal")
res_old_vs_ad <- lfcShrink(data.analisis, contrast = c("condition", "AD", "Old"), 
                           res = res_old_vs_ad, type = "normal")

summary(res_old_vs_ad)
summary(res_old_vs_young)
summary(res_ad_vs_young)


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

# Highlight first 5 PCs in one color, rest in dark grey
pc_data <- data.frame(prop_variance, pc = 1:length(prop_variance))
pc_data$highlight <- ifelse(pc_data$pc <= 5, "PC1-PC5 (~52% of variance)", "Remaining PCs")

ggplot(pc_data, aes(x = pc, y = prop_variance, fill = highlight)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("PC1-PC5 (~52% of variance)" = "darkorange", "Remaining PCs" = "darkturquoise")) +
  scale_y_continuous(limits = c(0, 0.25)) +
  scale_x_continuous(breaks = seq(0, length(prop_variance), by = 5)) +
  labs(x = "Principal components", y = "Proportion of variance explained", fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")

# Sample-to-sample correlation heatmap
pheatmap(vsd_cor_data)

# Save all three pairwise DESeq2 results for downstream scripts
saveRDS(res_old_vs_ad, "data/res_data.rds")
saveRDS(res_old_vs_young, "data/res_old_vs_young.rds")
saveRDS(res_ad_vs_young, "data/res_ad_vs_young.rds")

# Save data for use in downstream scripts
saveRDS(col.data, "data/col_data.rds")
