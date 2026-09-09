from pathlib import Path
from classes.Molecule import Molecule, Residu, Atom
from vedo import Points, Sphere, Plotter, Axes
import argparse
from typing import Dict, Any
import numpy as np


def getArgs() -> Dict[str, Any]:
    """
    Parses command-line arguments and returns them as a dictionary.
    
    Returns:
        dict: A dictionary containing the parsed command-line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--pdb_file", required=True, type=str)
    parser.add_argument("-n", "--n_sphere_points", required=True, type=int)

    args = parser.parse_args()
    input = Path(args.pdb_file)

    if not input.exists():
        raise ValueError("Input file does not exist")

    return {
        "input": input,
        "n_sphere_points": int(args.n_sphere_points),
    }


if __name__ == "__main__":

    args = getArgs()
    file = args.get("input")
    n_points = args.get("n_sphere_points", 92)

    molecule = Molecule.from_pdb_file(file, n_points)   
    molecule.compute_accessible_points()

    accessible_points = []
    for atom in molecule.atoms:
        accessible_points.extend(atom.accessible_points_list)
    
    pts = Points(accessible_points, r=10, c="red")
    spheres = [Sphere(pos=atom.coords, r=atom.radius + Atom.vdw_radius.get("water"))
                .alpha(0.15)
                .c("grey") for atom in molecule.atoms]

    coords = np.array([atom.coords for atom in molecule.atoms])

    xmin, ymin, zmin = coords.min(axis=0)
    xmax, ymax, zmax = coords.max(axis=0)

    axes = Axes(xrange=(xmin, xmax), yrange=(ymin, ymax), zrange=(zmin, zmax))
    plotter = Plotter(axes = axes, bg = "white")
    plotter.show(spheres, pts, axes = 1, viewup = "z", interactive = True)