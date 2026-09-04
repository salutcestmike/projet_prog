from classes.Atom import Atom
import numpy as np


class Residu:

    def __init__(self, resnum, resname, chain, atoms = None):
        self.resnum = resnum
        self.resname = resname
        self.chain = chain
        self.atoms = atoms

    def __str__(self):
        return f"{self.resnum} / {self.resname} / {self.chain} / {len(self.atoms)}"

    def get_accessible_surface(self, n_points = None, chain = None):
        if chain is not None and self.chain != chain:
            return 0
        return np.sum([atom.get_accessible_surface(n_points) for atom in self.atoms])

    def get_max_surface(self, chain = None):
        if chain is not None and self.chain != chain:
            return 0
        return np.sum([4 * np.pi * ((atom.radius + Atom.vdw_radius.get("water")) ** 2) for atom in self.atoms])