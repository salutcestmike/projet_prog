#!/bin/bash

# Créer l'environnement
conda create -y -n projet_prog python=3.11

# Installer les dépendances
conda install -y numpy tqdm scipy vedo