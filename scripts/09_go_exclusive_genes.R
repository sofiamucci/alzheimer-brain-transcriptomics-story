library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# == AD versus Old exclusive genes ========
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


# == OLd versus Young exclusive genes ========
# Load exclusive Old vs Young gene set
exclusive_old_young_df <- readRDS("data/genes_exclusive_old_vs_young.rds")

# Strip version suffix
gene_ids_clean_oy <- sub("\\..*", "", exclusive_old_young_df$gene_id)

# GO over-representation analysis
go_exclusive_old_young <- enrichGO(gene = gene_ids_clean_oy,
                                   OrgDb = org.Hs.eg.db,
                                   keyType = "ENSEMBL",
                                   ont = "ALL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.1,
                                   readable = TRUE)

as.data.frame(go_exclusive_old_young)


# == AD versus Young exclusive genes ========

# Load exclusive Old vs Young gene set
exclusive_ad_young_df <- readRDS("data/genes_exclusive_ad_vs_young.rds")

# Strip version suffix
gene_ids_clean_ay <- sub("\\..*", "", exclusive_ad_young_df$gene_id)

# GO over-representation analysis
go_exclusive_ad_young <- enrichGO(gene = gene_ids_clean_ay,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENSEMBL",
                                  ont = "ALL",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.1,
                                  readable = TRUE)

as.data.frame(go_exclusive_ad_young)
