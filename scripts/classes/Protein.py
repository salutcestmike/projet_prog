from classes.Atom import Atom, eucl_dist
from tqdm import tqdm
import numpy as np


class Protein:

    def __init__(self, prot_name, atoms, n_points):
        self.prot_name = prot_name
        self.n_points = n_points
        self.atoms = self.generate_atoms(atoms)
        self.accessible_points = 0
        self.neighbor_table = self.generate_neighbor_table()
        self.inaccessible_points = 0
        self.accessible_surface = 0

    def generate_atoms(self, atoms):
        return [Atom(int(atom[6:11]),
                     atom[12:16].strip(), 
                     atom[17:20].strip(),
                     int(atom[23:26]), 
                     [float(atom[30:38]), float(atom[39:46]), float(atom[47:54])]) for atom in atoms]

    def __str__(self):
        chaine = ""
        for atom in self.atoms:
            chaine += str(atom) + "\n"
        return chaine

    def generate_neighbor_table(self):
        table = np.zeros((len(self.atoms), len(self.atoms))).astype(bool)
        for i in tqdm(range(len(self.atoms))):
            for j in range(i + 1, len(self.atoms)):
                table[i][j] = eucl_dist(self.atoms[i].coords, self.atoms[j].coords) < (2 * Atom.vdw_radius.get("water") + self.atoms[i].radius + self.atoms[j].radius)
                table[j][i] = table[i][j]
        self.neighbor_table = table
        return table
    
    def count_inaccessible_points(self, n_points = None):
        if n_points is None:
            n_points = self.n_points
        self.inaccessible_points = 0

        for i in tqdm(range(len(self.atoms))):
            is_covered_table = np.zeros(n_points).astype(bool)
            sphere_points = self.atoms[i].get_points_v2(n_points)
            for j in range(len(self.atoms)):
                if (self.neighbor_table[i, j]):
                    for k in range(len(sphere_points)):
                        if eucl_dist(sphere_points[k], self.atoms[j].coords) < self.atoms[j].radius + Atom.vdw_radius.get("water"):
                            is_covered_table[k] = True
                            continue
            self.atoms[i].inaccessible_points = np.count_nonzero(is_covered_table)
            for j in range(len(sphere_points)):
                if is_covered_table[j]:
                    self.atoms[i].inaccessible_points_list.append(sphere_points[j])
                else:
                    self.atoms[i].accessible_points_list.append(sphere_points[j])
            self.inaccessible_points += np.count_nonzero(is_covered_table)

        return self.inaccessible_points
    
    def get_accessible_surface(self, n_points = None):
        if n_points is None:
            n_points = self.n_points
            
        self.accessible_surface = 0

        for atom in self.atoms:
            self.accessible_surface += (n_points - atom.inaccessible_points) * 4 * np.pi * (atom.radius ** 2) / n_points
        return self.accessible_surface

    def get_max_surface(self):
        max_surface = 0
        for atom in self.atoms:
            max_surface += 4 * np.pi * (atom.radius ** 2) 
        return max_surface