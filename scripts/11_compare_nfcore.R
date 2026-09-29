# 11_compare_nfcore.R
# Compare nf-core/rnaseq output against the manual pipeline output.
#
# nf-core/rnaseq command used to generate the comparison data:
#   nextflow run nf-core/rnaseq -profile docker -c custom.config \
#     --input samplesheet.csv --outdir /data/results/rnaseq_alzheimer \
#     --pseudo_aligner salmon \
#     --fasta .../GRCh38.primary_assembly.genome.fa.gz (GENCODE v47) \
#     --gtf .../gencode.v47.primary_assembly.annotation.gtf.gz (GENCODE v47)
#   custom.config: params { skip_alignment = true }  (salmon pseudo-alignment only)


##-----------------------------------------------------------------------------
# Step 1: align nf-core/rnaseq output with the manual pipeline output by gene_id
##-----------------------------------------------------------------------------

library(dplyr)

## Load nf-core/rnaseq output
# Adjust the path below to wherever this file lives on your machine
nfcore_counts <- read.delim("data/nfcore_comparison/salmon.merged.gene_counts.tsv")

## Load manual pipeline output (counts + TPM in one object)
txi <- readRDS("data/txi_full.rds")

manual_counts <- as.data.frame(txi$counts)
manual_tpm    <- as.data.frame(txi$abundance)

# gene_id is a rowname here, not a column — pull it out as a column
manual_counts$gene_id <- rownames(manual_counts)
manual_tpm$gene_id    <- rownames(manual_tpm)

## Strip the GENCODE version suffix from gene_id on both sides
# e.g. "ENSG00000000003.16" -> "ENSG00000000003"
# (same approach already used in scripts/08_venn_overlap_and_heatmap.R)
nfcore_counts$gene_id_clean <- sub("\\..*", "", nfcore_counts$gene_id)
manual_counts$gene_id_clean <- sub("\\..*", "", manual_counts$gene_id)
manual_tpm$gene_id_clean    <- sub("\\..*", "", manual_tpm$gene_id)

## Sanity check: how many genes on each side, any duplicates after stripping version?
cat("nf-core genes (raw):", nrow(nfcore_counts), "\n")
cat("nf-core genes (unique after stripping version):", n_distinct(nfcore_counts$gene_id_clean), "\n")
cat("nf-core duplicated gene_id_clean:", sum(duplicated(nfcore_counts$gene_id_clean)), "\n\n")

cat("manual genes (raw):", nrow(manual_counts), "\n")
cat("manual genes (unique after stripping version):", n_distinct(manual_counts$gene_id_clean), "\n")
cat("manual duplicated gene_id_clean:", sum(duplicated(manual_counts$gene_id_clean)), "\n")

## Join counts on the version-stripped gene_id
# suffix disambiguates overlapping sample-name columns from each source
joined_counts <- inner_join(
  manual_counts, nfcore_counts,
  by = "gene_id_clean",
  suffix = c("_manual", "_nfcore")
)

cat("\nGenes after inner join (counts):", nrow(joined_counts), "\n")

##-----------------------------------------------------------------------------
# Step 2: correlate counts between manual pipeline and nf-core, per sample
##-----------------------------------------------------------------------------

library(ggplot2)

# Sample IDs (SRR accessions) present in both matrices — same as the column
# names before the join added the _manual/_nfcore suffixes
sample_ids <- setdiff(colnames(manual_counts), c("gene_id", "gene_id_clean"))

## Correlation per sample, on log2(count + 1) scale
cor_results <- lapply(sample_ids, function(s) {
  manual_col <- joined_counts[[paste0(s, "_manual")]]
  nfcore_col <- joined_counts[[paste0(s, "_nfcore")]]
  
  log_manual <- log2(manual_col + 1)
  log_nfcore <- log2(nfcore_col + 1)
  
  data.frame(
    sample   = s,
    pearson  = cor(log_manual, log_nfcore, method = "pearson"),
    spearman = cor(log_manual, log_nfcore, method = "spearman")
  )
})

cor_summary <- do.call(rbind, cor_results)
cor_summary

## Overall summary across all samples
cat("Pearson  — mean:", mean(cor_summary$pearson),  " range:", range(cor_summary$pearson), "\n")
cat("Spearman — mean:", mean(cor_summary$spearman), " range:", range(cor_summary$spearman), "\n")

## Scatter plot for one representative sample (the first one)
# Swap sample_ids[1] for a specific SRR if you want to inspect a particular one
example_sample <- sample_ids[1]

plot_data <- data.frame(
  manual = log2(joined_counts[[paste0(example_sample, "_manual")]] + 1),
  nfcore = log2(joined_counts[[paste0(example_sample, "_nfcore")]] + 1)
)

ggplot(plot_data, aes(x = manual, y = nfcore)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Manual pipeline: log2(count + 1)",
    y = "nf-core/rnaseq: log2(count + 1)",
    title = paste("Counts comparison —", example_sample)
  )

ggsave("figures/nfcore_counts_comparison.png", width = 6, height = 5, dpi = 150)

#Step 2b: identify the genes driving the biggest divergence between
# manual and nf-core counts, for the same representative sample

# Residual from the y = x diagonal, in log2 space
plot_data$gene_id_clean <- joined_counts$gene_id_clean
plot_data$gene_name     <- joined_counts$gene_name   # nf-core already provides this
plot_data$diff          <- plot_data$manual - plot_data$nfcore

# Sort by absolute difference, biggest divergence first
top_divergent <- plot_data[order(-abs(plot_data$diff)), ]

# Show the top 20, with both raw values for context
head(top_divergent[, c("gene_id_clean", "gene_name", "manual", "nfcore", "diff")], 20)

# Step 2c: re-check correlation restricted to protein-coding genes only,
# using the same biotype logic already used in 07_volcano_plot.R

gene_biotype <- read.csv("data/gene_biotype_map.csv")
protein_coding_ids <- sub("\\..*", "", gene_biotype$gene_id[gene_biotype$biotype == "protein_coding"])

# Filter the joined counts table to protein-coding genes only, then re-run
# the same per-sample correlation loop as before
joined_pc <- joined_counts %>% filter(gene_id_clean %in% protein_coding_ids)

cat("Genes before filtering:", nrow(joined_counts), "\n")
cat("Genes after protein-coding filter:", nrow(joined_pc), "\n")

cor_results_pc <- lapply(sample_ids, function(s) {
  log_manual <- log2(joined_pc[[paste0(s, "_manual")]] + 1)
  log_nfcore <- log2(joined_pc[[paste0(s, "_nfcore")]] + 1)
  
  data.frame(
    sample   = s,
    pearson  = cor(log_manual, log_nfcore, method = "pearson"),
    spearman = cor(log_manual, log_nfcore, method = "spearman")
  )
})

cor_summary_pc <- do.call(rbind, cor_results_pc)

cat("Protein-coding only — Pearson  mean:", mean(cor_summary_pc$pearson),  "\n")
cat("Protein-coding only — Spearman mean:", mean(cor_summary_pc$spearman), "\n")


##-----------------------------------------------------------------------------
# TPM comparison — mirrors the counts comparison above
##-----------------------------------------------------------------------------

# Load nf-core TPM (already copied to data/nfcore_comparison/, same as counts)
nfcore_tpm <- read.delim("data/nfcore_comparison/salmon.merged.gene_tpm.tsv")
nfcore_tpm$gene_id_clean <- sub("\\..*", "", nfcore_tpm$gene_id)

# Join on version-stripped gene_id
joined_tpm <- inner_join(
  manual_tpm, nfcore_tpm,
  by = "gene_id_clean",
  suffix = c("_manual", "_nfcore")
)
cat("Genes after inner join (TPM):", nrow(joined_tpm), "\n")

# Restrict to protein-coding, same filter as counts
joined_tpm_pc <- joined_tpm %>% filter(gene_id_clean %in% protein_coding_ids)
cat("Genes after protein-coding filter (TPM):", nrow(joined_tpm_pc), "\n")

# Correlation per sample, log2(TPM + 1) scale — same reasoning as counts:
# log-transform so highly-expressed genes don't dominate the correlation
cor_results_tpm <- lapply(sample_ids, function(s) {
  log_manual <- log2(joined_tpm_pc[[paste0(s, "_manual")]] + 1)
  log_nfcore <- log2(joined_tpm_pc[[paste0(s, "_nfcore")]] + 1)
  data.frame(
    sample   = s,
    pearson  = cor(log_manual, log_nfcore, method = "pearson"),
    spearman = cor(log_manual, log_nfcore, method = "spearman")
  )
})
cor_summary_tpm <- do.call(rbind, cor_results_tpm)

cat("TPM, protein-coding — Pearson  mean:", mean(cor_summary_tpm$pearson),  "\n")
cat("TPM, protein-coding — Spearman mean:", mean(cor_summary_tpm$spearman), "\n")

