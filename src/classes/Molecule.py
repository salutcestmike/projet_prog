from classes.Residu import Residu, Atom
from tqdm import tqdm
import numpy as np
from numba import njit
from typing import List, Optional, Tuple, Union
from pathlib import Path


@njit
def generate_neighbor_table_numba(coords: np.ndarray, radii: np.ndarray, water_radius: float) -> np.ndarray:
    """
    Generates a neighbor table for the atoms in the molecule. The neighbor table is a boolean matrix where each entry (i, j) indicates whether atom i and atom j are neighbors based on their coordinates and radii.

    Args:
        coords (np.ndarray): An array of shape (n, 3) containing the coordinates of the atoms.
        radii (np.ndarray): An array of shape (n,) containing the radii of the atoms.
        water_radius (float): The van der Waals radius of water molecules.

    Returns:
        np.ndarray: A boolean matrix representing the neighbor relationships between atoms in the molecule.
    """
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

    def __init__(self, name : str, n_points : int, residues : Optional[List[Residu]] = None):
        self.name = name
        self.n_points = n_points
        self.residues = residues
        self.atoms = [atom for residue in self.residues for atom in residue.atoms]
        self.neighbor_table = self.generate_neighbor_table()

    @classmethod
    def from_pdb_file(cls, filename : Union[str, Path], n_points : int) -> 'Molecule':
        """
        Creates a Molecule instance from a PDB file by parsing the file to extract atom and residue information.
        
        Args:
            filename (Union[str, Path]): The path to the PDB file.
            n_points (int): The number of points to generate on each atom's surface.

        Returns:
            'Molecule': An instance of the Molecule class.
        """
        atoms_lines = []
        with open(filename, "r") as f:
            for line in f:
                if line.startswith(("ATOM")) and (line[77].strip() != "H"):
                    atoms_lines.append(line.strip())

        residues = [[]]
        resnum = 1
        resname = None
        chain = None
        for atom in atoms_lines:
            if int(atom[22:26]) != resnum:
                residues[len(residues) - 1] = Residu(resnum, resname, atom[21], residues[len(residues) - 1])
                residues.append([])
                resnum += 1
            if chain is not None and chain != atom[21]:
                resnum = 1
            resname = atom[17:20].strip()
            chain = atom[21]
            residues[len(residues) - 1].append(Atom(int(atom[6:11]),
                                                atom[12:16].strip(), 
                                                resname,
                                                [float(atom[30:38]), float(atom[39:46]), float(atom[47:54])]))
        residues[len(residues) - 1] = Residu(resnum, resname, atom[21], residues[len(residues) - 1])

        return cls(filename.stem, n_points, residues)

    def __str__(self) -> str:
        chaine = ""
        for residue in self.residues:
            chaine += str(residue) + "\n"
        return chaine

    def generate_neighbor_table(self) -> np.ndarray:
        """
        Generates a neighbor table for the atoms in the molecule. The neighbor table is a boolean matrix where each entry (i, j) indicates whether atom i and atom j are neighbors based on their coordinates and radii.
        
        Returns:
            numpy.ndarray: A boolean matrix representing the neighbor relationships between atoms in the molecule."""
        coords = np.asarray([atom.coords for atom in self.atoms], dtype=np.float64)
        radii = np.asarray([atom.radius for atom in self.atoms], dtype=np.float64)
        return generate_neighbor_table_numba(coords, radii, Atom.vdw_radius["water"])

    def compute_accessible_points(self, n_points: Optional[int] = None) -> Tuple[int, int]:
        """
        Computes the number of accessible and inaccessible points for each atom in the molecule based on their coordinates, radii, and the number of points generated on their surfaces.
        
        Returns:
            Tuple[int, int]: A tuple containing the count of accessible and inaccessible points.
        """
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
                if (self.neighbor_table[i][j]):
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

    def get_accessible_surface(self, n_points : Optional[int] = None, chain : Optional[str] = None) -> float:
        """
        Computes the accessible surface area of the molecule by summing the accessible surface areas of its residues.

        Args:
            n_points (int): The total number of points generated on each atom's surface.
            chain (str): The chain to which the residues belong.
        
        Returns:
            float: The accessible surface area of the molecule.
        """
        if n_points is None:
            n_points = self.n_points

        for atom in self.atoms:
            if (atom.accessible_points_list is None) or (atom.inaccessible_points_list is None):
                self.compute_accessible_points(n_points)
                break
        
        return sum([residue.get_accessible_surface(n_points, chain) for residue in self.residues])

    def get_max_surface(self, chain : Optional[str] = None) -> float:
        """
        Computes the maximum surface area of the molecule by summing the maximum surface areas of its residues.
            
        Args:
            chain (str): The chain to which the residues belong.
        
        Returns:
            float: The maximum surface area of the molecule.
        """
        return sum([residue.get_max_surface(chain) for residue in self.residues])