#!/bin/bash
# 01_download_and_quantify.sh
# Downloads, pseudo-aligns (salmon), and cleans up each sample sequentially
# to stay within available disk space.

cd ~/Documents/GitHub/alzheimer-brain-transcriptomics-story/data/fastq

samples=("SRR6145415" "SRR6145416" "SRR6145417" "SRR6145418" "SRR6145419"  
         "SRR6145420" "SRR6145421" "SRR6145422" "SRR6145423" "SRR6145424"  
         "SRR6145425" "SRR6145426" "SRR6145427" "SRR6145428" "SRR6145429"  
         "SRR6145430" "SRR6145431" "SRR6145432" "SRR6145433" "SRR6145434"  
         "SRR6145435" "SRR6145436" "SRR6145437" "SRR6145438" "SRR6145439" 
         "SRR6145440" "SRR6145441" "SRR6145442" "SRR6145443" "SRR6145444")

for srr in "${samples[@]}"; do
  echo "=== Processing $srr ==="
  
  prefetch "$srr"
  fasterq-dump "$srr"
  
  salmon quant -i ~/bioinformatics_references/GENCODE_v47/salmon_index \
    -l A \
    -r "${srr}.fastq" \
    -p 4 \
    --validateMappings \
    -o "${srr}_quant"
  
  rm "${srr}.fastq"
  rm -rf "$srr"
  
  echo "=== Done with $srr ==="
done

echo "All samples processed!"