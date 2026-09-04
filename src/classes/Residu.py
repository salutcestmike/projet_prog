from classes.Atom import Atom

class Residu:

    def __init__(self, resnum, resname, chain, atom_list = None):
        self.resnum = resnum
        self.resname = resname
        self.chain = chain
        self.atom_list = atom_list

    def __str__(self):
        return f"{self.resnum} / {self.resname} / {self.chain} / {len(self.atom_list)}"
    