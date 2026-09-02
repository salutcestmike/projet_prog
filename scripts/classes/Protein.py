import classes.Atom 

class Protein:

    def __init__(self, prot_name, atoms):
        self.prot_name = prot_name
        self.atoms = self.generate_atoms(atoms)
        
    def generate_atoms(self, atoms):
        return [classes.Atom (atom[6:11],
                              atom[12:16], 
                              atom[17:20],
                              atom[23:26], 
                              [float(atom[30:38]), float(atom[39:46]), float(atom[47:54])]) for atom in atoms]

