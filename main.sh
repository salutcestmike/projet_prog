#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate projet_prog
python -u ./scripts/main.py
conda deactivate
