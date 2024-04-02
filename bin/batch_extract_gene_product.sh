#!/bin/bash

# Path to the folder containing GenBank files
input_folder="/mnt/scratch/runs/Pangraph_test/ALL_input/Indiv_input_fasta/Prokka"

# Iterate over each GenBank file in the folder
for file in "$input_folder"/*.gbk; do
    # Extract file name without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"

    # Run the Python script for each GenBank file
    python extract_sequences.py "$file" -p "IS1380 family transposase ISEcp1" "Beta-lactamase CTX-M-1"
done
