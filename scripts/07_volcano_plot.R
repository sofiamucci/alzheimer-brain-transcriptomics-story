library(DESeq2)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(Biostrings)

# Build gene_id → biotype mapping from GENCODE fasta headers
fasta_headers <- names(readDNAStringSet("~/bioinformatics_references/GENCODE_v47/gencode.v47.transcripts.fa.gz"))
header_parts <- strsplit(fasta_headers, "\\|")

gene_biotype <- data.frame(
  gene_id = sapply(header_parts, `[`, 2),
  biotype = sapply(header_parts, function(x) x[length(x)])
)
gene_biotype <- unique(gene_biotype)

protein_coding_ids <- sub("\\..*", "", gene_biotype$gene_id[gene_biotype$biotype == "protein_coding"])

# Load all three contrasts (generated in 05_qc.R, padj<0.05, shrunk)
res_ad_vs_old <- readRDS("data/res_data.rds")
res_old_vs_young <- readRDS("data/res_old_vs_young.rds")
res_ad_vs_young <- readRDS("data/res_ad_vs_young.rds")

# == AD vs Old  ==============================================
# ============================================================
res_table_ad_old <- as.data.frame(res_ad_vs_old)
res_table_ad_old <- tibble::rownames_to_column(res_table_ad_old, "gene_id")

res_table_ad_old <- res_table_ad_old %>% 
  mutate(Expression = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Unchanged"
  ))

# Filter to protein-coding genes only (applied post-hoc — see METHODS.md)
res_table_ad_old <- res_table_ad_old %>%
  filter(sub("\\..*", "", gene_id) %in% protein_coding_ids)

table(res_table_ad_old$Expression)

# Top genes for labeling
top_genes_ad_old <- res_table_ad_old %>%
  filter(Expression != "Unchanged") %>%
  arrange(padj) %>%
  head(30)
top_genes_ad_old$gene_id_clean <- sub("\\..*", "", top_genes_ad_old$gene_id)
top_genes_ad_old$gene_name <- mapIds(org.Hs.eg.db, keys = top_genes_ad_old$gene_id_clean, 
                                     keytype = "ENSEMBL", column = "SYMBOL")
top_genes_ad_old$label <- ifelse(is.na(top_genes_ad_old$gene_name), 
                                 top_genes_ad_old$gene_id_clean, 
                                 top_genes_ad_old$gene_name)

res_table_ad_old <- left_join(res_table_ad_old, top_genes_ad_old[, c("gene_id", "label")], by = "gene_id")

# Volcano plot
p_ad_old <- ggplot(res_table_ad_old, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col = Expression)) + 
  scale_color_manual(values = c("Down-regulated" = "darkturquoise", 
                                "Unchanged" = "grey25", 
                                "Up-regulated" = "darkorange")) +
  geom_vline(xintercept = c(-1, 1), colour = "grey", linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), colour = "grey", linetype = "dotted", linewidth = 0.5) +
  ggtitle("AD vs Old") +
  theme_bw() +
  theme(legend.position = "top")

p_ad_old + 
  ggrepel::geom_text_repel(
    aes(label = label), 
    size = 3, 
    max.overlaps = 30,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey50"
  )

# == Old vs Young  ===========================================
# ============================================================
res_table_old_young <- as.data.frame(res_old_vs_young)
res_table_old_young <- tibble::rownames_to_column(res_table_old_young, "gene_id")

res_table_old_young <- res_table_old_young %>% 
  mutate(Expression = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Unchanged"
  ))

# Filter to protein-coding genes only
res_table_old_young <- res_table_old_young %>%
  filter(sub("\\..*", "", gene_id) %in% protein_coding_ids)

table(res_table_old_young$Expression)

# Top genes for labeling 
top_genes_old_young <- res_table_old_young %>%
  filter(Expression != "Unchanged") %>%
  arrange(padj) %>%
  head(30)
top_genes_old_young$gene_id_clean <- sub("\\..*", "", top_genes_old_young$gene_id)
top_genes_old_young$gene_name <- mapIds(org.Hs.eg.db, keys = top_genes_old_young$gene_id_clean, 
                                        keytype = "ENSEMBL", column = "SYMBOL")
top_genes_old_young$label <- ifelse(is.na(top_genes_old_young$gene_name), 
                                    top_genes_old_young$gene_id_clean, 
                                    top_genes_old_young$gene_name)

res_table_old_young <- left_join(res_table_old_young, top_genes_old_young[, c("gene_id", "label")], by = "gene_id")

# Volcano plot
p_old_young <- ggplot(res_table_old_young, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col = Expression)) + 
  scale_color_manual(values = c("Down-regulated" = "darkturquoise", 
                                "Unchanged" = "grey25", 
                                "Up-regulated" = "darkorange")) +
  geom_vline(xintercept = c(-1, 1), colour = "grey", linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), colour = "grey", linetype = "dotted", linewidth = 0.5) +
  ggtitle("Old vs Young") +
  theme_bw() +
  theme(legend.position = "top")

p_old_young + 
  ggrepel::geom_text_repel(
    aes(label = label), 
    size = 3, 
    max.overlaps = 30,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey50"
  )

# == AD vs Young  ============================================
# ============================================================
res_table_ad_young <- as.data.frame(res_ad_vs_young)
res_table_ad_young <- tibble::rownames_to_column(res_table_ad_young, "gene_id")

res_table_ad_young <- res_table_ad_young %>% 
  mutate(Expression = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Unchanged"
  ))

# Filter to protein-coding genes only
res_table_ad_young <- res_table_ad_young %>%
  filter(sub("\\..*", "", gene_id) %in% protein_coding_ids)

table(res_table_ad_young$Expression)

# Top genes for labeling
top_genes_ad_young <- res_table_ad_young %>%
  filter(Expression != "Unchanged") %>%
  arrange(padj) %>%
  head(30)
top_genes_ad_young$gene_id_clean <- sub("\\..*", "", top_genes_ad_young$gene_id)
top_genes_ad_young$gene_name <- mapIds(org.Hs.eg.db, keys = top_genes_ad_young$gene_id_clean, 
                                       keytype = "ENSEMBL", column = "SYMBOL")
top_genes_ad_young$label <- ifelse(is.na(top_genes_ad_young$gene_name), 
                                   top_genes_ad_young$gene_id_clean, 
                                   top_genes_ad_young$gene_name)

res_table_ad_young <- left_join(res_table_ad_young, top_genes_ad_young[, c("gene_id", "label")], by = "gene_id")

# Volcano plot
p_ad_young <- ggplot(res_table_ad_young, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col = Expression)) + 
  scale_color_manual(values = c("Down-regulated" = "darkturquoise", 
                                "Unchanged" = "grey25", 
                                "Up-regulated" = "darkorange")) +
  geom_vline(xintercept = c(-1, 1), colour = "grey", linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), colour = "grey", linetype = "dotted", linewidth = 0.5) +
  ggtitle("AD vs Young") +
  theme_bw() +
  theme(legend.position = "top")

p_ad_young + 
  ggrepel::geom_text_repel(
    aes(label = label), 
    size = 3, 
    max.overlaps = 30,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey50"
  )

#== Save data ===============================================
saveRDS(res_table_ad_old, "data/res_table_ad_old_final.rds")
saveRDS(res_table_old_young, "data/res_table_old_young_final.rds")
saveRDS(res_table_ad_young, "data/res_table_ad_young_final.rds")
