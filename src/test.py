import numpy as np
import freesasa

def surface(atom1, atom2, radius1, radius2):
    dist = np.sqrt((atom1[0] - atom2[0])** 2 + 
                    (atom1[1] - atom2[1])** 2 + 
                    (atom1[2] - atom2[2])** 2)
    h1 = radius1 - (((dist ** 2) + (radius1 ** 2) - (radius2 ** 2)) / (2 * dist))
    h2 = radius2 - (((dist ** 2) + (radius2 ** 2) - (radius1 ** 2)) / (2 * dist))
    surface = 4 * np.pi * (radius1 ** 2) + 4 * np.pi * (radius2 ** 2) - 2 * np.pi * radius1 * h1 - 2 * np.pi * radius2 * h2  
    return surface

atom1 = [0, 0, 0]
atom2 = [2, 0, 0]
radius1 = 2
radius2 = 2

print(surface(atom1, atom2, radius1, radius2))


structure = freesasa.Structure("data/pdb/lysozyme_253L.pdb")
structure = freesasa.Structure("data/pdb/benzene.pdb")
result = freesasa.calc(structure)
area_classes = freesasa.classifyResults(result, structure)

print("Total : %.2f A2" % result.totalArea())
for key in area_classes:
    print(key, ": %.2f A2" % area_classes[key])
