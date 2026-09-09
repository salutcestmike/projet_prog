from pathlib import Path


with open(Path("data/pdb_ids.txt"), "r") as f:
    ids = [line.strip() for line in f.readlines()]
print(ids)
print(len(ids))
# results = Path("output/results")
# naccess_results = Path("output/naccess_results")

# resultats = []
# for sous_dossier in results.iterdir():
#     if not sous_dossier.is_dir():
#         continue

#     identifiant = sous_dossier.name

#     # Calcul de la surface
#     surface = calculer_surface(sous_dossier)

#     resultats.append({
#         "identifiant": identifiant,
#         "surface": surface
#     })

# df = pd.DataFrame(resultats)

# print(df)