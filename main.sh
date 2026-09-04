#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate projet_prog

python3 -u ./src/main.py \
    -i data/pdb/3I40.pdb \
    -o output/results \
    -n 92 

conda deactivate
