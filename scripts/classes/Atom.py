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
        self.aromatic = resname in ["HIS", "TRP", "TYR", "PHE"]
        self.radius = self.get_radius()
        self.inaccessible_points = 0
        self.accessible_points_list = []
        self.inaccessible_points_list = []

    def __str__(self):
        return f"{self.radius} / {self.atom_name} / {self.atom_num} / {self.resname} / {self.resnum} / coords({self.coords[0]}, {self.coords[1]}, {self.coords[2]})"

    def get_radius(self):
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

    def get_points_v1(self, n_points = 92):
        radius = self.radius + Atom.vdw_radius.get("water")

        points = []
        
        golden_angle = np.pi * (3 - np.sqrt(5))
    
        for i in range(n_points):
    
            z = (1 - 2 * (i + 0.5) / n_points)
            r = np.sqrt(1 - z*z) * radius 
    
            theta = golden_angle * i
    
            x = r * np.cos(theta)
            y = r * np.sin(theta)
    
            points.append([x + self.coords[0], y + self.coords[1], z * radius + self.coords[2]])
    
        return np.array(points)

    def get_points_v2(self, n_points = 92):
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

def eucl_dist(coords1, coords2):
    return np.sqrt((coords1[0] - coords2[0]) ** 2 + (coords1[1] - coords2[1]) ** 2 + (coords1[2] - coords2[2]) ** 2)