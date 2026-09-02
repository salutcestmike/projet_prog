import numpy as np
import json


class Atom:

    with open('config/vdw_radius.json') as f:
        vdw_radius = json.load(f)

    def __init__(self, atom_num, atom_name, resname, resnum, coords):
        self.atom_num = atom_num
        self.atom_name = atom_name
        self.resname = resname
        self.resnum = resnum
        self.coords = coords
        self.aromatic = resname in ["HIS", "TRP", "TYR", "PHE"] # KEEP ?
        self.radius = self.get_radius()

    def __str__(self):
        return f"{self.radius} / {self.atom_name} / {self.atom_num} / {self.resname} / {self.resnum} / coords({self.coords[0]}, {self.coords[1]}, {self.coords[2]})"

    def get_radius(self):
        print(self.atom_name[0])
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
                if (self.resname == "PHE" and self.atom_name in ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]) or 
                   (self.resname == "PHE" and self.atom_name in ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]) or 
                   (self.resname == "PHE" and self.atom_name in ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]) or 
                   (self.resname == "PHE" and self.atom_name in ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]): 

                
            return 2
            # if self.aromatic:
        return None
             