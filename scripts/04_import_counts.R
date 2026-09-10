# Load libraries
library(tximport)
library(dplyr)

# Load sample metadata
sample_info <- read.csv("data/sample_metadata_full.csv")

# Build paths to each sample's quant.sf file
quant_files <- file.path("data/fastq", paste0(sample_info$external_id, "_quant"), "quant.sf")
names(quant_files) <- sample_info$external_id

# Confirm all files exist before importing
all(file.exists(quant_files))

# Build transcript-to-gene mapping from GENCODE fasta headers
library(Biostrings)

fasta_headers <- names(readDNAStringSet("~/bioinformatics_references/GENCODE_v47/gencode.v47.transcripts.fa.gz"))

tx2gene <- data.frame(
  TXNAME = sapply(strsplit(fasta_headers, "\\|"), `[`, 1),
  GENEID = sapply(strsplit(fasta_headers, "\\|"), `[`, 2)
)

head(tx2gene)

# Import salmon quantifications, aggregated to gene level
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# Extract raw counts matrix
counts <- txi$counts
dim(counts)
head(counts[, 1:5])

# Save files
saveRDS(txi, "data/txi_full.rds")
write.csv(counts, "data/counts_full.csv", row.names = TRUE)


