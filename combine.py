import numpy as np
from pymatgen.core import Structure
from sys import argv
from latticeMismatch import *

struc_1 = Structure.from_file(argv[1])
struc_2 = Structure.from_file(argv[2])


print()
print()

l1 = struc_1.lattice.matrix
l2 = struc_2.lattice.matrix

print("l1: ")
print(l1)
print()


print("l2: ")
print(l2)
print()

# Selecting the 2D matrix of the strucs
l1_2D = [[l1[0,0], l1[0,1]], [l1[1,0], l1[1,1]]]
l2_2D = [[l2[0,0], l2[0,1]], [l2[1,0], l2[1,1]]]

print("l1_2D: ")
print(l1_2D)
print()

print("l2_2D: ")
print(l2_2D)
print()

"""
Calculates superlattices

returns [bestM, bestN, bestR, smallestErr]

struc_1 * lattice2 * rotation
struc_2 * lattice1

"""

transformations = latticeTransformations(l1_2D, l2_2D, 20, 0.01)

m = transformations[0]
n = transformations[1]
r = transformations[2]

m_3D = [[m[0][0], m[0][1], 0], [m[1][0], m[1][1], 0], [0, 0, 1]]
n_3D = [[n[0][0], n[0][1], 0], [n[1][0], n[1][1], 0], [0, 0, 1]]

print(m_3D)

struc_1.make_supercell(scaling_matrix=n_3D)
struc_1.rotate_sites(theta=np.radians(r))

struc_2.make_supercell(scaling_matrix=m_3D)

print("supercell 1: ")
print(struc_1.lattice)
print()


print("supercell 2: ")
print(struc_2.lattice)
print()