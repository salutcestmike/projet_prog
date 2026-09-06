from classes.Atom import Atom
import numpy as np
from typing import List, Optional

class Residu:

    def __init__(self, resnum : int, resname : str, atoms : List[Atom], chain : Optional[str] = None,):
        self.resnum = resnum
        self.resname = resname
        self.chain = chain
        self.atoms = atoms

    def __str__(self) -> str:
        return f"{self.resnum} / {self.resname} / {self.chain} / {len(self.atoms)}"

    def get_accessible_surface(self, n_points : Optional[int] = None, chain : Optional[str] = None) -> float:
        """
        Computes the accessible surface area of the residue by summing the accessible surface areas of its atoms.
        
        Args:
            n_points (Optional[int]): The total number of points generated on each atom's surface.
            chain (Optional[str]): The chain to which the residue belongs.

        Returns:
            float: The accessible surface area of the residue.
        """
        # If a specific chain is provided and it does not match the residue's chain
        if chain is not None and self.chain != chain:
            return 0
        return np.sum([atom.get_accessible_surface(n_points) for atom in self.atoms])

    def get_max_surface(self, chain : Optional[str] = None) -> float:
        """
        Computes the maximum surface area of the residue by summing the maximum surface areas of its atoms.
        
        Args:
            chain (Optional[str]): The chain to which the residue belongs.
            
        Returns:
            float: The maximum surface area of the residue.
        """
        # If a specific chain is provided and it does not match the residue's chain
        if chain is not None and self.chain != chain:
            return 0
        return np.sum([atom.get_max_surface() for atom in self.atoms])

    def get_surface_results(self, n_points : Optional[int] = None) -> List[dict]:
        """
        Generates a list of dictionaries containing the surface results for each atom in the residue.
        
        Args:
            n_points (Optional[int]): The total number of points generated on each atom's surface.
            
        Returns:
            List[dict]: A list of dictionaries containing the atom number, atom name, and surfaces for each atom the residue.
        """
        results = []
        for atom in self.atoms:
            accessible_surface = atom.get_accessible_surface(n_points)
            max_surface = atom.get_max_surface()
            results.append({
                "atom_num": atom.atom_num,
                "atom_name": atom.atom_name,
                "accessible_surface": accessible_surface,
                "max_surface": max_surface,
                "accessible_percentage": (accessible_surface / max_surface) * 100,
            })
        return results