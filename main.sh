#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate projet_prog

python3 -u ./scripts/main.py \
    -i data/pdb/lysozyme_253L.pdb \
    -o output/surfaces \
    -w 3 \
    -n 92 \
    -d 

conda deactivate
