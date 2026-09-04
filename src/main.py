from pathlib import Path
from classes.Molecule import Molecule, Atom
from vedo import Points, Sphere, Plotter, Axes
import argparse
from datetime import datetime



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

    output_folder.mkdir(parents=True, exist_ok=True)

    return {
        "input": input,
        "output_folder": output_folder,
        "num_workers": int(args.num_workers),
        "n_sphere_points": int(args.n_sphere_points),
        "display_molecule": args.display_molecule if input.is_file() else False
    }



if __name__ == "__main__":

    args = getArgs()

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    output_folder = Path(f"{args.get("output_folder")}/{timestamp}")
    output_folder.mkdir(parents=True, exist_ok=True)

    files_to_process = []
    if args.get("input").is_file():
        files_to_process.append(args.get("input"))
    else:
        files_to_process = list(args.get("input").glob("*.pdb"))
    print(f"Files to process:")
    [print(f"{file}\n") for file in files_to_process]

    display = args.get("display_molecule", False)
    n_points = args.get("n_sphere_points", 92)
    for file in files_to_process:
        atoms = []
        atoms_type = []

        with open(file, "r") as f:
            for line in f:
                if line.startswith(("ATOM")) and (line[77].strip() != "H"):
                    atoms.append(line.strip())
                    atoms_type.append(line[12:16].strip())
    
        molecule = Molecule(file.stem, atoms, n_points)
        accessible_points_count, inaccessible_points_count = molecule.compute_accessible_points(n_points)

        print(f"Number of atoms: {len(molecule.atoms)}")
        print(f"Number of accessible points: {accessible_points_count}")
        print(f"Number of inaccessible points: {inaccessible_points_count}")
        print(f"{round(100 * accessible_points_count / (len(molecule.atoms) * n_points), 2)}% of points are accessible")
        print(f"Accessible surface: {molecule.get_accessible_surface()} Å2")
        print(f"Max surface: {round(molecule.get_max_surface(), 2)} Å2")

        with open(f"{output_folder}/{file.stem}.txt", 'w') as f:
            for atom in molecule.atoms:
                f.write(f"{atom.atom_num} / {atom.atom_name} / {100 * len(atom.accessible_points_list) / n_points}%\n")

        if display:
            display = False

            accessible_points = []
            for atom in molecule.atoms:
                accessible_points.extend(atom.accessible_points_list)
            
            pts = Points(
                accessible_points,
                r=10,
                c="red"
            )
        
            spheres = [Sphere(
                    pos=atom.coords,
                    r=atom.radius + Atom.vdw_radius.get("water")
                ).alpha(0.15).c("grey") for atom in molecule.atoms]
            
        
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