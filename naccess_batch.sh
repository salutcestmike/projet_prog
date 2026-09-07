#!/bin/bash

PDB_DIR="data/pdb"
OUTPUT_DIR="./output/naccess_results"

mkdir -p "$OUTPUT_DIR"

for pdb in ./"$PDB_DIR"/*.pdb; do

    filename=$(basename "$pdb" .pdb)

    echo "Processing $filename"

    mkdir -p "$OUTPUT_DIR/$filename"

    (
        cd "$OUTPUT_DIR/$filename" || exit 1
        echo "$pdb"
        naccess "../../../$PDB_DIR/$filename.pdb"

    )

done