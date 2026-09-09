#!/bin/bash
# 00_setup_environment.sh
# Sets up the conda environment and downloads reference files needed 
# for pseudoalignment with salmon.

# Create conda environment with salmon, sra-tools, and fastqc
conda create -n rnaseq -c bioconda -c conda-forge salmon sra-tools fastqc -y

# Activate environment
conda activate rnaseq

# Download human transcriptome reference (GENCODE v47)
mkdir -p reference
cd reference
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz

# Build salmon index
salmon index -t gencode.v47.transcripts.fa.gz -i salmon_index -k 31