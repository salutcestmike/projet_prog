import Atom

class Residu:
    count = 0

    def __init__(self, resname, atoms):
        self.resname = resname
        self.aromatic = True if resname in ["HIS", "TRP", "TYR", tyrosine,]
        self.atoms = atoms