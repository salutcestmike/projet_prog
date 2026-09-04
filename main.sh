#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate projet_prog

python3 -u ./src/main.py \
    -i data/pdb/1PMA.pdb \
    -o output/results \
    -w 8 \
    -n 92 

conda deactivate
