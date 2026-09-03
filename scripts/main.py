from pathlib import Path
from classes.Protein import Protein, Atom
from vedo import Points, Sphere, Plotter, Axes


if __name__ == "__main__":
    pdb_file_test1 = Path("data/pdb/test.pdb")
    pdb_file_test2 = Path("data/pdb/benzene.pdb")
    pdb_file1 = Path("data/pdb/insulin_3I40.pdb")
    pdb_file2 = Path("data/pdb/lysozyme_253L.pdb")
    pdb_file3 = Path("data/pdb/lysozyme_1GWD.pdb")

    atoms = []
    atoms_type = []

    with open(pdb_file_test2, "r") as f:
        for line in f:
            if line.startswith(("ATOM")):
                atoms.append(line.strip())
                atoms_type.append(line[12:16].strip())

    [print(a) for a in atoms]
    prot = Protein("insulin", atoms)
    print(prot)
    print(set(atoms_type))
    print(prot.neighbor_table)

    inaccessible_points = prot.count_inaccessible_points()
    accessible_surface = prot.get_accessible_surface()
    print(f"Number of atoms: {len(prot.atoms)}")
    print(f"Number of inaccessible points: {inaccessible_points}")
    print(f"{100 * inaccessible_points / (len(prot.atoms) * 92)}% of points are inaccessible")
    print(accessible_surface)

    points = []
    for atom in prot.atoms:
        points.extend(atom.accessible_points)

    pts = Points(
        points,
        r=10,
        c="red"
    )

    spheres = [Sphere(
        pos=atom.coords,
        r=atom.radius + Atom.vdw_radius.get("water")
    ).alpha(0.15).c("black" if (atom.atom_name[0] == "C") 
                            else ("blue" if (atom.atom_name[0] == "N") 
                                  else ("red" if (atom.atom_name[0] == "O") else "yellow"))) for atom in prot.atoms]

    spheres = [Sphere(
            pos=atom.coords,
            r=atom.radius + Atom.vdw_radius.get("water")
        ).alpha(0.15).c("grey") for atom in prot.atoms]
    

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
        pts,
        axes = 1,
        viewup="z",
        interactive=True
    )