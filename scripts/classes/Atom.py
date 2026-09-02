import numpy as np

class Atom:

    def __init__(self, atom_num, atom_name, resname, resnum, coords):
        self.atom_num = atom_num
        self.atom_name = atom_name
        self.resname = resname
        self.resnum = resnum
        self.coords = coords
        self.aromatic = resname in ["HIS", "TRP", "TYR", "PHE"]