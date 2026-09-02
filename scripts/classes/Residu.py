import classes.Atom 

class Residu:
    count = 0

    def __init__(self, resname, resnum, atoms):
        self.resname = resname
        self.resnum = resnum
        self.aromatic = resname in ["HIS", "TRP", "TYR", "PHE"]
        self.atoms = self.generate_atoms(atoms)

    def generate_atoms(self, atoms):
        return [classes.Atom (atom[6:11],
                              atom[12:16], 
                              [float(atom[30:38]), float(atom[39:46]), float(atom[47:54])]) for atom in atoms]

