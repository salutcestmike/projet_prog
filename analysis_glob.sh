#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate projet_prog

python3 -u ./src/analysis_glob.py  

conda deactivate
