# Data

This directory does not contain raw sequencing data or fastq files. Raw 
reads are downloaded and processed programmatically via 
`scripts/02_download_and_quantify.sh`, which retrieves data from GEO/SRA 
and pseudo-aligns with salmon.

## Source
- **GEO accession**: GSE104704
- **SRA study**: SRP119561
- **Samples**: all 30 samples (8 Young, 10 Old, 12 AD)

## Contents
- `sample_metadata_full.csv` — sample IDs and condition labels
- `normalized_data.rds/csv` — DESeq2-normalized counts
- `res_data.rds`, `res_old_vs_young.rds`, `res_ad_vs_young.rds` — DESeq2 
  results per contrast
- `res_table_sig_*.rds` — significant genes per contrast (protein-coding, 
  padj<0.05, |log2FC|>0.5)
- `genes_exclusive_*.rds/csv` — genes exclusive to each contrast
- `txi_full.rds` — tximport object with gene-level quantifications

## Reproducing the pipeline
Run scripts `00` through `10` in `scripts/` in numbered order. See 
[METHODS.md](../METHODS.md) for full methodological detail.