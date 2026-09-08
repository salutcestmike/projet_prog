from pathlib import Path
from classes.Molecule import Molecule, Residu, Atom
from vedo import Points, Sphere, Plotter, Axes
import argparse
from datetime import datetime
from typing import Dict, Any
from tqdm import tqdm



def getArgs() -> Dict[str, Any]:
    """
    Parses command-line arguments and returns them as a dictionary.
    
    Returns:
        dict: A dictionary containing the parsed command-line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--pdb_folder_or_file", required=True, type=str)
    parser.add_argument("-o", "--output_folder", required=True, type=str)
    parser.add_argument("-n", "--n_sphere_points", required=True, type=int)

    args = parser.parse_args()

    input = Path(args.pdb_folder_or_file)
    output_folder = Path(args.output_folder)

    if not input.exists():
        raise ValueError("Input folder/file does not exist")

    output_folder.mkdir(parents=True, exist_ok=True)

    return {
        "input": input,
        "output_folder": output_folder,
        "n_sphere_points": int(args.n_sphere_points),
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
    [print(f"{file}") for file in files_to_process]
    print()

    n_points = args.get("n_sphere_points", 92)

    for file in tqdm(files_to_process, desc="Files analyzed"):
        molecule = Molecule.from_pdb_file(file, n_points)   

        accessible_points_count, inaccessible_points_count = molecule.compute_accessible_points()

        print(f"Molecule: {molecule.name}")
        print(f"Number of residues: {len(molecule.residues)}")
        print(f"Number of atoms: {len(molecule.atoms)}")
        print(f"Number of accessible points: {accessible_points_count}")
        print(f"Number of inaccessible points: {inaccessible_points_count}")

        chains = list(set([residue.chain for residue in molecule.residues]))
        if len(chains) > 1:
            chains.sort()
            for chain in chains:
                accessible_surface = molecule.get_accessible_surface(chain)
                max_accessible_surface = molecule.get_max_surface(chain)
                print(f"Chain {chain} accessible surface: {round(accessible_surface, 2)} Å^2")
                print(f"Chain {chain} accessible percentage: {round(100 * accessible_surface / max_accessible_surface, 2)} %")

        print(f"Total accessible surface: {round(molecule.get_accessible_surface(), 2)} Å^2")
        print(f"Total accessible percentage: {round(100 * accessible_points_count / (len(molecule.atoms) * n_points), 2)} %")
        print()

        # Export results
        molecule.export_surface_results(output_folder)
        