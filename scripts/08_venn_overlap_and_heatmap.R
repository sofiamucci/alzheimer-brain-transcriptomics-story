library(dplyr)
library(ggVennDiagram)
library(org.Hs.eg.db)
library(VennDiagram)


# Load final (protein-coding filtered) results from 07_volcano_plot.R
res_table_ad_old <- readRDS("data/res_table_ad_old_final.rds")
res_table_old_young <- readRDS("data/res_table_old_young_final.rds")
res_table_ad_young <- readRDS("data/res_table_ad_young_final.rds")

# Extract and save significant genes for each contrast
res_table_sig_ad_old <- res_table_ad_old %>% filter(Expression != "Unchanged")
res_table_sig_old_young <- res_table_old_young %>% filter(Expression != "Unchanged")
res_table_sig_ad_young <- res_table_ad_young %>% filter(Expression != "Unchanged")

saveRDS(res_table_sig_ad_old, "data/res_table_sig_ad_old.rds")
saveRDS(res_table_sig_old_young, "data/res_table_sig_old_young.rds")
saveRDS(res_table_sig_ad_young, "data/res_table_sig_ad_young.rds")

nrow(res_table_sig_ad_old)
nrow(res_table_sig_old_young)
nrow(res_table_sig_ad_young)

# Venn diagram of overlap across the three comparisons
gene_lists <- list(
  "AD vs Old" = res_table_sig_ad_old$gene_id,
  "Old vs Young" = res_table_sig_old_young$gene_id,
  "AD vs Young" = res_table_sig_ad_young$gene_id
)

venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL,
  fill = c("darkorange", "mediumpurple", "darkturquoise"),
  alpha = 0.6,
  col = "black",
  lwd = 1,
  cex = 1.2,
  cat.cex = 1.1,
  cat.col = "black",
  print.mode = c("raw", "percent"),
  sigdigits = 2,
  main = "Overlap of significant genes across comparisons"
)

grid::grid.newpage()
grid::grid.draw(venn.plot)

# Extract gene lists per comparison
genes_ad_old <- res_table_sig_ad_old$gene_id
genes_old_young <- res_table_sig_old_young$gene_id
genes_ad_young <- res_table_sig_ad_young$gene_id

# Genes shared between each pair
shared_ad_old_and_old_young <- intersect(genes_ad_old, genes_old_young)
shared_ad_old_and_ad_young <- intersect(genes_ad_old, genes_ad_young)
shared_old_young_and_ad_young <- intersect(genes_old_young, genes_ad_young)

# --- Genes shared across all three comparisons ---
shared_all_three <- Reduce(intersect, list(genes_ad_old, genes_old_young, genes_ad_young))

# --- Function to map Ensembl IDs to gene symbols ---
map_to_symbol <- function(gene_ids) {
  clean_ids <- sub("\\..*", "", gene_ids)
  mapIds(org.Hs.eg.db, keys = clean_ids, keytype = "ENSEMBL", column = "SYMBOL")
}

# --- Print counts and gene names ---
length(shared_ad_old_and_old_young)
map_to_symbol(shared_ad_old_and_old_young)

length(shared_ad_old_and_ad_young)
map_to_symbol(shared_ad_old_and_ad_young)

length(shared_old_young_and_ad_young)
map_to_symbol(shared_old_young_and_ad_young)

length(shared_all_three)
map_to_symbol(shared_all_three)

# --- Genes exclusive to AD vs Old (not in either of the other two contrasts) ---
exclusive_ad_old <- setdiff(genes_ad_old, union(genes_old_young, genes_ad_young))

length(exclusive_ad_old)
map_to_symbol(exclusive_ad_old)

exclusive_ad_old_df <- data.frame(
  gene_id = exclusive_ad_old,
  gene_symbol = map_to_symbol(exclusive_ad_old)
)

saveRDS(exclusive_ad_old_df, "data/genes_exclusive_ad_vs_old.rds")
write.csv(exclusive_ad_old_df, "data/genes_exclusive_ad_vs_old.csv", row.names = FALSE)

exclusive_ad_old_df

# --- Genes exclusive to each contrast (not shared with the other two) ---
exclusive_old_young <- setdiff(genes_old_young, union(genes_ad_old, genes_ad_young))
length(exclusive_old_young)
map_to_symbol(exclusive_old_young)

exclusive_ad_young <- setdiff(genes_ad_young, union(genes_ad_old, genes_old_young))
length(exclusive_ad_young)
map_to_symbol(exclusive_ad_young)

exclusive_old_young_df <- data.frame(
  gene_id = exclusive_old_young,
  gene_symbol = map_to_symbol(exclusive_old_young)
)

exclusive_ad_young_df <- data.frame(
  gene_id = exclusive_ad_young,
  gene_symbol = map_to_symbol(exclusive_ad_young)
)

saveRDS(exclusive_old_young_df, "data/genes_exclusive_old_vs_young.rds")
write.csv(exclusive_old_young_df, "data/genes_exclusive_old_vs_young.csv", row.names = FALSE)

saveRDS(exclusive_ad_young_df, "data/genes_exclusive_ad_vs_young.rds")
write.csv(exclusive_ad_young_df, "data/genes_exclusive_ad_vs_young.csv", row.names = FALSE)


# Combined heatmap: significant genes across all three contrasts
# ============================================================
library(pheatmap)

normalized_data <- readRDS("data/normalized_data.rds")
col.data <- readRDS("data/col_data.rds")

# --- Union of significant genes across all three contrasts ---
all_sig_genes <- unique(c(
  res_table_sig_ad_old$gene_id,
  res_table_sig_old_young$gene_id,
  res_table_sig_ad_young$gene_id
))

length(all_sig_genes)

# --- Subset normalized counts to these genes ---
heatmap_data <- normalized_data[all_sig_genes, ]

# --- Annotation for column (sample) colors by condition ---
annotation_col <- data.frame(condition = col.data$condition)
rownames(annotation_col) <- rownames(col.data)

annotation_colors <- list(
  condition = c("Young" = "darkturquoise", "Old" = "grey40", "AD" = "darkorange")
)

# --- Heatmap ---
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cutree_rows = 3,
         main = "Significant genes across all three contrasts (n=30 samples)")
