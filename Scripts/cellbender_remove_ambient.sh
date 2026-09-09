#!/bin/bash

#1. Get the directory where THIS script is saved
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#2. Point to your input and output directories dynamically
INPUT_DIR="${SCRIPT_DIR}/raw_count_tables/cellbender_input"
OUTPUT_DIR="${SCRIPT_DIR}/clean_adata"

#3. Create the output folder if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Starting CellBender Pipeline"
echo "Input folder: $INPUT_DIR"

#4. Loop through every .h5ad file in the folder
for h5ad_file in "$INPUT_DIR"/*.h5ad; do
    
    #Check if files actually exist to prevent loop errors
    [ -e "$h5ad_file" ] || continue
    
    #Extract just the sample name (e.g., "Sample_A" from "/path/Sample_A.h5ad")
    FILENAME=$(basename "$h5ad_file")
    SAMPLE_NAME="${FILENAME%.h5ad}"
    
    echo "Processing Sample: $SAMPLE_NAME"
    
    #5. Run the core CellBender command
    #NOTE: If you are running locally on your Mac CPU, remove the '--cuda' flag!
    cellbender remove-background \
        --input "$h5ad_file" \
        --output "${OUTPUT_DIR}/${SAMPLE_NAME}_cellbender.h5" \
        --total-droplets-included 50000 \
        --fpr 0.01 \
        --epochs 150 \
        --cuda

done

echo "All samples successfully processed by CellBender"

