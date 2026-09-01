from pathlib import Path


pdb_file1 = Path("data/pdb/insulin_3I40.pdb")
pdb_file2 = Path("data/pdb/lysozyme_253L.pdb")

atom_names = []

with open(pdb_file1, "r") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16].strip()
            atom_names.append(atom_name)

print(set(atom_names))
