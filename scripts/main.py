from pathlib import Path
from classes.Molecule import Molecule, Atom
from vedo import Points, Sphere, Plotter, Axes
import argparse



def getArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--pdb_folder_or_file", required=True, type=str)
    parser.add_argument("-o", "--output_folder", required=True, type=str)
    parser.add_argument("-w", "--num_workers", required=True, type=int)
    parser.add_argument("-n", "--n_sphere_points", required=True, type=int)
    parser.add_argument("-d", "--display_molecule", action="store_true")

    args = parser.parse_args()

    input = Path(args.pdb_folder_or_file)
    output_folder = Path(args.output_folder)

    if not input.exists():
        raise ValueError("Input folder/file does not exist")

    output_folder.parent.mkdir(parents=True, exist_ok=True)

    return {
        "input": input,
        "output_file": output_folder,
        "num_workers": int(args.num_workers),
        "n_sphere_points": int(args.n_sphere_points),
        "display_molecule": args.display_molecule if input.is_file() else False
    }


if __name__ == "__main__":

    args = getArgs()

    files_to_process = []
    if args.get("input").is_file():
        files_to_process.append(args.get("input"))
    else:
        files_to_process = list(args.get("input").glob("*.pdb"))
    print(f"Files to process:")
    [print(f"{file}\n") for file in files_to_process]


    n_points = 92
    pdb_file_test1 = Path("data/pdb/test.pdb")
    pdb_file_test2 = Path("data/pdb/benzene.pdb")
    pdb_file1 = Path("data/pdb/insulin_3I40.pdb")
    pdb_file2 = Path("data/pdb/lysozyme_253L.pdb")
    pdb_file3 = Path("data/pdb/lysozyme_1GWD.pdb")
    pdb_file4 = Path("data/pdb/1MH1.pdb")

    atoms = []
    atoms_type = []

    with open(pdb_file2, "r") as f:
        for line in f:
            if line.startswith(("ATOM")):
                atoms.append(line.strip())
                atoms_type.append(line[12:16].strip())

    [print(a) for a in atoms]
    prot = Molecule("insulin", atoms, n_points)
    print(prot)
    print(set(atoms_type))
    print(prot.neighbor_table)

    inaccessible_points = prot.count_inaccessible_points()
    accessible_surface = prot.get_accessible_surface()
    print(f"Number of atoms: {len(prot.atoms)}")
    print(f"Number of inaccessible points: {inaccessible_points}")
    print(f"{100 * inaccessible_points / (len(prot.atoms) * n_points)}% of points are inaccessible")
    print(accessible_surface)
    print(prot.get_max_surface())

    with open("res.txt", 'w') as f:
        for atom in prot.atoms:
            f.write(f"{atom.atom_num} / {atom.atom_name} / {100 * len(atom.accessible_points_list) / n_points}%\n")
    
    points = []
    for atom in prot.atoms:
        points.extend(atom.accessible_points_list)

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
