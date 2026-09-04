from classes.Atom import Atom, eucl_dist
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool, cpu_count

class Molecule:

    def __init__(self, prot_name, atoms, n_points):
        self.prot_name = prot_name
        self.n_points = n_points
        self.atoms = self.generate_atoms(atoms)
        self.neighbor_table = self.generate_neighbor_table()

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
        return table

    def compute_accessible_points(self, n_points = None):
        if n_points is None:
            n_points = self.n_points

        accessible_points_count = 0
        inaccessible_points_count = 0

        for i in tqdm(range(len(self.atoms))):
            self.atoms[i].accessible_points_list = []
            self.atoms[i].inaccessible_points_list = []
            is_covered_table = np.zeros(n_points).astype(bool)
            sphere_points = np.array(self.atoms[i].get_points(n_points))

            for j in range(len(self.atoms)):
                if (self.neighbor_table[i, j]):
                    for k in range(len(sphere_points)):
                        if eucl_dist(sphere_points[k], self.atoms[j].coords) < self.atoms[j].radius + Atom.vdw_radius.get("water"):
                            is_covered_table[k] = True
                            continue

            for j in range(len(sphere_points)):
                if is_covered_table[j]:
                    self.atoms[i].inaccessible_points_list.append(sphere_points[j])
                else:
                    self.atoms[i].accessible_points_list.append(sphere_points[j])

            accessible_points_count += len(self.atoms[i].accessible_points_list)
            inaccessible_points_count += len(self.atoms[i].inaccessible_points_list)

        return accessible_points_count, inaccessible_points_count

    def get_accessible_surface(self, n_points = None):
        if n_points is None:
            n_points = self.n_points

        for atom in self.atoms:
            if (atom.accessible_points_list is None) or (atom.inaccessible_points_list is None):
                self.compute_accessible_points(n_points)
                break
        
        accessible_surface = 0

        for atom in self.atoms:
            accessible_surface += atom.get_accessible_surface(n_points)
        return accessible_surface

    def get_max_surface(self):
        max_surface = 0
        for atom in self.atoms:
            max_surface += 4 * np.pi * ((atom.radius + Atom.vdw_radius.get("water")) ** 2) 
        return max_surface