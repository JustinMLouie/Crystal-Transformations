import numpy as np
from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import RotationTransformation
from sys import argv
from latticeMismatch import *

struc_1 = Structure.from_file(argv[1])
struc_2 = Structure.from_file(argv[2])

print()
print()

l1 = struc_1.lattice.matrix
l2 = struc_2.lattice.matrix

# Selecting the 2D matrix of the strucs
l1_2D = l1[:2,:2]
l2_2D = l2[:2,:2]

"""
Calculates superlattices

returns [bestM, bestN, bestR, smallestErr]

struc_1 * lattice2 * rotation
struc_2 * lattice1
"""

transformations = latticeTransformations(l1_2D, l2_2D, 100, 0.01)

m = transformations[0]
n = transformations[1]
r = transformations[2]

m_3D = [[m[0][0], m[0][1], 0],
		 [m[1][0], m[1][1], 0],
		  [0, 0, 1]]

n_3D = [[n[0][0], n[0][1], 0],
		 [n[1][0], n[1][1], 0], 
		  [0, 0, 1]]

struc_1.make_supercell(scaling_matrix=n_3D)
struc_1 = RotationTransformation([0,0,1], r, False).apply_transformation(struc_1)

struc_2.make_supercell(scaling_matrix=m_3D)

struc_1_coords = struc_1.cart_coords
struc_2_coords = struc_2.cart_coords

"""
Repositions struc_2 directly above struc_1


Finds the highest Z in struc_1
Finds the lowest Z in struc_2
Moves struc_2 to be arg[v] above struc_1
"""

maxZ1 = 0
minZ2 = 1000

for particle in struc_1_coords:
	if (particle[2] >= maxZ1):
		maxZ1 = particle[2]

for particle in struc_2_coords:
	if (particle[2] <= minZ2):
		minZ2 = particle[2]

zAdjust = maxZ1 - minZ2 + float(argv[3])

for particle in struc_2_coords:
	particle[2] += zAdjust 

# struc_1.to(filename='POSCAR_{}'.format('superlattice'))



