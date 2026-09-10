from classes.Atom import Atom
import numpy as np
from typing import List, Optional
import json


class Residu:

    with open('config/max_asa.json') as f:
        max_asa = json.load(f)

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
        Computes the maximum surface area of the residue defined by the config/max_asa.json file.
        
        Args:
            chain (Optional[str]): The chain to which the residue belongs.
            
        Returns:
            float: The maximum surface area of the residue.
        """
        # If a specific chain is provided and it does not match the residue's chain
        if chain is not None and self.chain != chain:
            return 0
        return Residu.max_asa.get(self.resname, 0)

    def get_surface_results(self, n_points : Optional[int] = None) -> dict:
        """
        Generates a dictionary containing the surface results for the residue, including its number, name, chain, and the surface results of its atoms.
        
        Args:
            n_points (Optional[int]): The total number of points generated on each atom's surface.
            
        Returns:
            dict: A dictionary containing the residue number, residue name, chain, and a list of
            surface results for each atom in the residue.
        """
        return {
            "resnum": self.resnum,
            "resname": self.resname,
            "chain": self.chain,
            "atoms": [atom.get_surface_results(n_points) for atom in self.atoms]
        }