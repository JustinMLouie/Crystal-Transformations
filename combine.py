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

print("l1: ")
print(l1)
print()


print("l2: ")
print(l2)
print()

# Selecting the 2D matrix of the strucs
l1_2D = l1[:2,:2]
l2_2D = l2[:2,:2]

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

print("m_3D: ")
print(m_3D)

print("n_3D: ")
print(n_3D)

struc_1.make_supercell(scaling_matrix=n_3D)
struc_1 = RotationTransformation([0,0,1], r, False).apply_transformation(struc_1)

struc_2.make_supercell(scaling_matrix=m_3D)

print("supercell 1: ")
print(struc_1.lattice)
print()

print("supercell 2: ")
print(struc_2.lattice)
print

"""
TODO: 

find the lowest z value of the top layer
find the highest z value of the bottom layer
stack them using these z values with gap of user input
struc with bigger area is base layer

"""

# struc_1.to(filename='POSCAR_{}'.format('superlattice'))



