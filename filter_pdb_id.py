from pathlib import Path
import shutil

dossier = Path("output/results/20260908103650")
dossier_naccess = Path("output/naccess_results")

ids = [d.name for d in dossier.iterdir() if d.is_dir()]

ids_naccess = []

for sous_dossier in dossier_naccess.iterdir():
    if sous_dossier.is_dir():
        fichier_rsa = sous_dossier / ".rsa"
        print(fichier_rsa)
        if fichier_rsa.exists() and fichier_rsa.stat().st_size > 0:
            ids_naccess.append(sous_dossier.name)

common_ids = set(ids) & set(ids_naccess)
print(common_ids)
print(len(common_ids))

for sous_dossier in dossier.iterdir():
    if sous_dossier.is_dir() and sous_dossier.name not in common_ids:
        shutil.rmtree(sous_dossier)

with open("pdb_ids.txt", "w") as f:
    for i in list(common_ids):  
        f.write(f"{i}\n")