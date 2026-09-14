library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Load exclusive AD vs Old gene set
exclusive_ad_old_df <- readRDS("data/genes_exclusive_ad_vs_old.rds")

# Strip version suffix
gene_ids_clean <- sub("\\..*", "", exclusive_ad_old_df$gene_id)

# GO over-representation analysis
go_exclusive_ad_old <- enrichGO(gene = gene_ids_clean,
                                OrgDb = org.Hs.eg.db,
                                keyType = "ENSEMBL",
                                ont = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.1,
                                readable = TRUE)

as.data.frame(go_exclusive_ad_old)

# Dotplot
dotplot(go_exclusive_ad_old, showCategory = 15) + 
  ggtitle("GO enrichment: genes exclusive to AD vs Old")


