from classes.Residu import Residu, Atom
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool, cpu_count
from numba import njit


@njit
def generate_neighbor_table_numba(coords, radii, water_radius):
    n = len(coords)
    table = np.zeros((n, n), dtype=np.bool)

    for i in range(n):
        for j in range(i + 1, n):
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            dz = coords[i][2] - coords[j][2]
            distance = dx*dx + dy*dy + dz*dz
            threshold = 2 * water_radius + radii[i] + radii[j]

            table[i][j] = distance < threshold * threshold
            table[j][i] = table[i][j]

    return table

class Molecule:

    def __init__(self, name, atoms, n_points):
        self.name = name
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
        coords = np.asarray([atom.coords for atom in self.atoms], dtype=np.float64)
        radii = np.asarray([atom.radius for atom in self.atoms], dtype=np.float64)
        return generate_neighbor_table_numba(coords, radii, Atom.vdw_radius["water"])

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
                    distances = np.linalg.norm(sphere_points - np.array(self.atoms[j].coords), axis=1)
                    is_covered_table = is_covered_table | (distances < self.atoms[j].radius + Atom.vdw_radius["water"])

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