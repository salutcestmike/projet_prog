import numpy as np

class Atom:

    def __init__(self, atom_num, atom_name, coords):
        self.atom_num = atom_num
        self.atom_name = atom_name
        self.coords = coords