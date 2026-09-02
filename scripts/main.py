from pathlib import Path
from classes.Protein import Protein
from vedo import Points, Sphere, Plotter, Axes


if __name__ == "__main__":
    pdb_file_test = Path("data/pdb/benzene.pdb")
    pdb_file1 = Path("data/pdb/insulin_3I40.pdb")
    pdb_file2 = Path("data/pdb/lysozyme_253L.pdb")

    atoms = []
    atoms_type = []

    with open(pdb_file_test, "r") as f:
        for line in f:
            if line.startswith(("ATOM")):
                atoms.append(line.strip())
                atoms_type.append(line[12:16].strip())

    [print(a) for a in atoms]
    prot = Protein("insulin", atoms)
    print(prot)
    print(set(atoms_type))
    print(prot.neighbor_table)

    print(prot.count_inaccessible_points())
    print(prot.get_accessible_surface())

    spheres = [Sphere(
        pos=atom.coords,
        r=atom.radius
    ).alpha(0.15).c("black" if (atom.atom_name[0] == "C") 
                            else ("blue" if (atom.atom_name[0] == "N") 
                                  else ("red" if (atom.atom_name[0] == "O") else "yellow"))) for atom in prot.atoms]

    axes = Axes(
        xtitle="X",
        ytitle="Y",
        ztitle="Z"
    )

    plotter = Plotter(
        axes=axes,
        bg="white"
    )

    plotter.show(
        spheres,
        axes,
        viewup="z",
        interactive=True
    )