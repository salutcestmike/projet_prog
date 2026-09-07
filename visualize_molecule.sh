#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate projet_prog

python3 -u ./src/visualize_molecule.py \
    -i data/pdb_test/1PMA.pdb \
    -n 92

conda deactivate
