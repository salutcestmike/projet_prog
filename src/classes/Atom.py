import numpy as np
import json
from typing import Tuple, Optional

class Atom:

    with open('config/vdw_radius.json') as f:
        vdw_radius = json.load(f)

    def __init__(self, atom_num : int, atom_name : str, resname : str, coords : Tuple[float, float, float]):
        self.atom_num = atom_num
        self.atom_name = atom_name
        self.coords = coords
        self.resname = resname
        self.aromatic = resname in ["HIS", "TRP", "TYR", "PHE"]
        self.radius = self.get_radius()
        self.accessible_points_list = None
        self.inaccessible_points_list = None

    def __str__(self) -> str:
        return f"{self.radius} / {self.atom_name} / {self.atom_num} / {self.resname} / {self.resnum} / coords({self.coords[0]}, {self.coords[1]}, {self.coords[2]})"

    def get_radius(self) -> float:
        """
        Returns the van der Waals radius of the atom based on its name and residue type, according to Shrake article.
        
        Returns:
            float: The van der Waals radius of the atom.
        """
        if self.atom_name[0] == "N":
            return Atom.vdw_radius.get("nitrogen")
        elif self.atom_name[0] == "S":
            return Atom.vdw_radius.get("sulfur")
        elif self.atom_name[0] == "O":
            return Atom.vdw_radius.get("oxygen")
        elif self.atom_name == "C":
            return Atom.vdw_radius.get("carbonyl and other carbon")
        elif self.atom_name[0] == "C":
            if not(self.aromatic):
                return Atom.vdw_radius.get("non-aromatic carbon")
            else:
                if ((self.resname in ["PHE", "TYR"] and self.atom_name in ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]) or 
                   (self.resname == "TRP" and self.atom_name in ["CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]) or 
                   (self.resname == "HIS" and self.atom_name in ["CG", "CD2", "CE1"])): 
                    return Atom.vdw_radius.get("aromatic carbon")
                else:
                    return Atom.vdw_radius.get("non-aromatic carbon")
        return None

    def get_points(self, n_points : Optional[int] = 92) -> np.ndarray:
        """
        Generates points on the surface of a sphere centered at the atom's coordinates according to Saff and Kuijlaars algorithm.
        
        Args:
            n_points (Optional[int]): The number of points to generate.
            
        Returns:
            numpy.ndarray: An array of 3D points on the surface of the sphere.
        """
        radius = self.radius + Atom.vdw_radius.get("water")

        thetas = []
        phis = []
        points = []

        for k in range(1, n_points + 1):
            hk = -1 + 2 * ((k - 1) / (n_points - 1))
            theta = np.arccos(hk)

            if ((k == 1) or (k == n_points)):
                phi = 0
            else:
                phi = (phis[k - 2] + (3.6 / np.sqrt(n_points)) * (1 / np.sqrt(1 - (hk ** 2)))) % (2*np.pi)

            thetas.append(theta)
            phis.append(phi)
            points.append([np.sin(theta) * np.cos(phi) * radius + self.coords[0],
                        np.sin(theta) * np.sin(phi) * radius + self.coords[1], 
                        np.cos(theta) * radius + self.coords[2]])
        return np.array(points)

    def get_accessible_surface(self, n_points : Optional[int] = None) -> float:
        """
        Calculates the accessible surface area of the atom based on the number of accessible points.
        
        Args:
            n_points (Optional[int]): The total number of points generated on the atom's surface.
            
        Returns:
            float: The accessible surface area of the atom.
        """
        # If n_points is not provided, it defaults to the sum of accessible and inaccessible points.
        if n_points is None:
            n_points = len(self.accessible_points_list) + len(self.inaccessible_points_list)
        return ((len(self.accessible_points_list) * 4 * np.pi * ((self.radius + Atom.vdw_radius.get("water")) ** 2)) / n_points)

    def get_max_surface(self) -> float:
        """
        Calculates the maximum surface area of the atom based on its radius and the van der Waals radius of water.
        
        Returns:
            float: The maximum surface area of the atom.
        """
        return 4 * np.pi * ((self.radius + Atom.vdw_radius.get("water")) ** 2)
