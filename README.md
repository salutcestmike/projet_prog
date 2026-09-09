# projet_prog

Outil Python d'analyse et de visualisation de structures protéiques (fichiers **PDB**) : calcul de l'accessibilité au solvant par échantillonnage de points sur des sphères atomiques, et visualisation 3D de la molécule avec mise en évidence des points non couverts.

## Sommaire

- [Description](#description)
- [Structure du projet](#structure-du-projet)
- [Prérequis](#prérequis)
- [Installation](#installation)
- [Utilisation](#utilisation)
  - [Lancer l'analyse](#lancer-lanalyse)
  - [Visualiser une molécule](#visualiser-une-molécule)
- [Exemple](#exemple)
- [Dépendances](#dépendances)
- [Auteur](#auteur)

## Description

Ce projet prend en entrée une structure protéique au format `.pdb` ou un dossier contenant les fichiers `.pdb` et :

1. place un ensemble de points répartis uniformément sur la sphère de chaque atome (nombre de points paramétrable via `-n`) ;
2. détermine, pour chaque point, s'il est couvert (enfoui) ou non couvert (accessible) en fonction des atomes voisins, ce qui permet d'estimer la surface accessible au solvant (SASA) ;
3. exporte les résultats dans un dossier de sortie ;
4. permet de visualiser la molécule en 3D avec les points de sphère non couverts mis en évidence.

## Structure du projet

```
projet_prog/
├── config/                 # Fichiers de configuration
├── data/                   # Données d'entrée (structures PDB, jeux de test)
├── doc/                    # Documentation complémentaire
├── src/                    # Code source Python
│   ├── main.py              # Point d'entrée : calcul de l'accessibilité
│   └── visualize_molecule.py # Visualisation 3D de la molécule
├── env_install.sh          # Création de l'environnement conda
├── main.sh                 # Lance le calcul sur un exemple
├── visualize_molecule.sh   # Lance la visualisation sur un exemple
└── README.md
```

## Prérequis

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) ou Anaconda installé
- Python 3.11
- Un accès à un affichage graphique si vous utilisez la visualisation 3D (vedo)

## Installation

Clonez le dépôt puis créez l'environnement conda dédié :

```bash
git clone https://github.com/salutcestmike/projet_prog.git
cd projet_prog
bash env_install.sh
```

Ce script crée un environnement conda nommé `projet_prog` (Python 3.11) et installe les dépendances : `numpy`, `tqdm`, `scipy`, `vedo`.

## Utilisation

### Lancer l'analyse

```bash
conda activate projet_prog

python3 -u ./src/main.py \
    -i data/pdb_test/3I40.pdb \
    -o output/results \
    -n 92

conda deactivate
```

Arguments :

| Option | Description |
|--------|-------------|
| `-i`   | Chemin du fichier PDB en entrée ou du dossier contenant les fichiers PDB à analyser |
| `-o`   | Dossier de sortie pour les résultats |
| `-n`   | Nombre de points échantillonnés par sphère atomique |

Vous pouvez aussi simplement exécuter le script fourni :

```bash
bash main.sh
```

### Visualiser une molécule

Pour afficher la structure en 3D avec les points de sphère non couverts :

```bash
conda activate projet_prog

python3 -u ./src/visualize_molecule.py \
    -i data/pdb_test/3I40.pdb \
    -n 92

conda deactivate
```

Ou directement :

```bash
bash visualize_molecule.sh
```

## Exemple

Un fichier de test est fourni dans `data/pdb_test/3I40.pdb` pour essayer rapidement le projet sans avoir à fournir vos propres données.

## Dépendances

- [numpy](https://numpy.org/)
- [tqdm](https://github.com/tqdm/tqdm)
- [scipy](https://scipy.org/)
- [vedo](https://vedo.embl.es/)

## Auteur

Projet développé par [salutcestmike](https://github.com/salutcestmike).
