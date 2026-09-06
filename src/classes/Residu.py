from classes.Atom import Atom
import numpy as np
from typing import List, Optional

class Residu:

    def __init__(self, resnum : int, resname : str, chain : str, atoms : Optional[List[Atom]] = None):
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
            n_points (int): The total number of points generated on each atom's surface.
            chain (str): The chain to which the residue belongs.
            
        Returns:
            float: The accessible surface area of the residue.
        """
        if chain is not None and self.chain != chain:
            return 0
        return np.sum([atom.get_accessible_surface(n_points) for atom in self.atoms])

    def get_max_surface(self, chain : Optional[str] = None) -> float:
        """
        Computes the maximum surface area of the residue by summing the maximum surface areas of its atoms.
        
        Args:
            chain (str): The chain to which the residue belongs.
            
        Returns:
            float: The maximum surface area of the residue.
        """
        if chain is not None and self.chain != chain:
            return 0
        return np.sum([4 * np.pi * ((atom.radius + Atom.vdw_radius.get("water")) ** 2) for atom in self.atoms])