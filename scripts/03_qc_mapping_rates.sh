#!/bin/bash
# 02_qc_mapping_rates.sh
# Checks mapping rates across all quantified samples to confirm the 
# pseudoalignment step ran successfully before proceeding to R analysis.

cd ~/Documents/GitHub/alzheimer-brain-transcriptomics-story/data/fastq

echo "=== Number of quant folders found ==="
ls -d SRR*_quant | wc -l

echo ""
echo "=== Mapping rate per sample ==="
for dir in SRR*_quant; do
  srr=$(basename "$dir" _quant)
  rate=$(grep -I "mapping rate" "$dir/logs/salmon_quant.log")
  echo "$srr: $rate"
done