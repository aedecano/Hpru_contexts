#!/bin/bash

# Path to the folder containing GenBank files
input_folder="/mnt/scratch/runs/Pangraph_test/ALL_input/Indiv_input_fasta/Prokka/"

# Iterate over each GenBank file in the folder
for file in "$input_folder"/*.gbk; do
    # Extract file name without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"

    # Run the Python script for each GenBank file
    python3 extract_startgene_endgene_collapse.py "$file" -s "IS1380 family transposase ISEcp1" -e "Beta-lactamase CTX-M-1"

    # Replace "hypothetical protein_1" with "spacer_proteein_1" in the output FASTA file
    sed -i 's/hypothetical protein_1/spacer_protein_1/g' "${filename_no_ext}_cut.fasta"
done
