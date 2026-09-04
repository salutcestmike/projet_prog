from classes.Residu import Residu, Atom
from tqdm import tqdm
import numpy as np
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

    def __init__(self, name, n_points, residues = None):
        self.name = name
        self.n_points = n_points
        self.residues = residues
        self.neighbor_table = self.generate_neighbor_table()

    def __str__(self):
        chaine = ""
        for residue in self.residues:
            chaine += str(residue) + "\n"
        return chaine

    def generate_neighbor_table(self):
        coords = np.asarray([atom.coords for residue in self.residues for atom in residue.atoms], dtype=np.float64)
        radii = np.asarray([atom.radius for residue in self.residues for atom in residue.atoms], dtype=np.float64)
        return generate_neighbor_table_numba(coords, radii, Atom.vdw_radius["water"])

    def compute_accessible_points(self, n_points = None):
        if n_points is None:
            n_points = self.n_points

        accessible_points_count = 0
        inaccessible_points_count = 0

        atoms = [atom for residue in self.residues for atom in residue.atoms]

        for i in tqdm(range(len(atoms))):
            atoms[i].accessible_points_list = []
            atoms[i].inaccessible_points_list = []
            is_covered_table = np.zeros(n_points).astype(bool)
            sphere_points = np.array(atoms[i].get_points(n_points))

            for j in range(len(atoms)):
                if (self.neighbor_table[i][j]):
                    distances = np.linalg.norm(sphere_points - np.array(atoms[j].coords), axis=1)
                    is_covered_table = is_covered_table | (distances < atoms[j].radius + Atom.vdw_radius["water"])

            for j in range(len(sphere_points)):
                if is_covered_table[j]:
                    atoms[i].inaccessible_points_list.append(sphere_points[j])
                else:
                    atoms[i].accessible_points_list.append(sphere_points[j])

            accessible_points_count += len(atoms[i].accessible_points_list)
            inaccessible_points_count += len(atoms[i].inaccessible_points_list)

        return accessible_points_count, inaccessible_points_count

    def get_accessible_surface(self, n_points = None, chain = None):
        if n_points is None:
            n_points = self.n_points

        for residue in self.residues:
            for atom in residue.atoms:
                if (atom.accessible_points_list is None) or (atom.inaccessible_points_list is None):
                    self.compute_accessible_points(n_points)
                    break
        
        return sum([residue.get_accessible_surface(n_points, chain) for residue in self.residues])

    def get_max_surface(self, chain = None):
        return sum([residue.get_max_surface(chain) for residue in self.residues])