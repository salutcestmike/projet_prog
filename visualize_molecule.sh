#!/bin/bash

# Script pour visualiser une protéine en 3D avec les points de sphère non couverts

source ~/miniconda3/etc/profile.d/conda.sh

conda activate projet_prog

python3 -u ./src/visualize_molecule.py \
    -i data/pdb_test/3I40.pdb \
    -n 92

conda deactivate
