library(DESeq2)
library(dplyr)
library(ggplot2)
library(pROC)
library(org.Hs.eg.db)

# Load data
normalized_data <- readRDS("data/normalized_data.rds")
col.data <- readRDS("data/col_data.rds")
res_table_sig_ad_old <- readRDS("data/res_table_sig_ad_old.rds")
exclusive_ad_old_df <- readRDS("data/genes_exclusive_ad_vs_old.rds")

# Candidates: genes from the immune/inflammatory GO categories 
# identified in 09_go_exclusive_genes.R
immune_genes <- c("MPO", "SELE", "ELANE", "PRTN3", "IL4I1", "KMO", "IL33", "CITED1", "EGR1")

candidates <- res_table_sig_ad_old %>%
  filter(sub("\\..*", "", gene_id) %in% sub("\\..*", "", exclusive_ad_old_df$gene_id))

candidates$gene_id_clean <- sub("\\..*", "", candidates$gene_id)
candidates$gene_symbol <- mapIds(org.Hs.eg.db, keys = candidates$gene_id_clean, 
                                 keytype = "ENSEMBL", column = "SYMBOL")

candidates <- candidates %>% filter(gene_symbol %in% immune_genes)

candidates[, c("gene_id", "gene_symbol", "log2FoldChange", "padj")]

# Prepare data for boxplots
candidate_counts <- normalized_data[candidates$gene_id, ]
rownames(candidate_counts) <- candidates$gene_symbol

candidate_long <- as.data.frame(candidate_counts) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "normalized_count")

candidate_long$condition <- col.data$condition[match(candidate_long$sample, rownames(col.data))]
candidate_long <- candidate_long %>% filter(condition %in% c("Old", "AD"))

# Boxplots
ggplot(candidate_long, aes(x = condition, y = normalized_count, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.7, position = position_jitter(width = 0.05, seed = 42)) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Old" = "darkturquoise", "AD" = "darkorange")) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(x = NULL, y = "Normalized counts", 
       title = "Immune/inflammatory candidate genes (AD vs Old, exclusive)")

# ROC/AUC per candidate gene
# Exploratory: genes were selected based on significance in this same 
# dataset (circular evaluation), so AUC values should not be interpreted 
# as validated biomarker performance — see note below.

roc_results <- list()
for (i in seq_len(nrow(candidates))) {
  gene_data <- candidate_long %>% filter(gene == candidates$gene_symbol[i])
  roc_obj <- roc(gene_data$condition, gene_data$normalized_count, quiet = TRUE)
  roc_results[[candidates$gene_symbol[i]]] <- data.frame(
    gene = candidates$gene_symbol[i], 
    auc = as.numeric(auc(roc_obj))
  )
}

roc_summary <- do.call(rbind, roc_results)
roc_summary <- roc_summary[order(-roc_summary$auc), ]
roc_summary






