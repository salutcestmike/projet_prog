from pathlib import Path
from classes.Protein import Protein


if __name__ == "__main__":
    pdb_file1 = Path("data/pdb/insulin_3I40.pdb")
    pdb_file2 = Path("data/pdb/lysozyme_253L.pdb")

    atoms = []
    atoms_type = []

    with open(pdb_file1, "r") as f:
        for line in f:
            if line.startswith(("ATOM")):
                atoms.append(line.strip())
                atoms_type.append(line[12:16].strip())

    [print(a) for a in atoms]
    prot = Protein("insulin", atoms)
    print(prot)
    print(set(atoms_type))