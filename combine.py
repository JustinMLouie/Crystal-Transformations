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

transformations = latticeTransformations(l1_2D, l2_2D, 20, 0.01)

m = transformations[0]
n = transformations[1]
r = transformations[2]

m_3D = [[m[0][0], m[0][1], 0],
		 [m[1][0], m[1][1], 0],
		  [0, 0, 1]]

n_3D = [[n[0][0], n[0][1], 0],
		 [n[1][0], n[1][1], 0], 
		  [0, 0, 1]]

struc_1.make_supercell(scaling_matrix=m_3D)
RotationTransformation([0,0,1], r, False).apply_transformation(struc_1)

struc_2.make_supercell(scaling_matrix=n_3D)

l1 = struc_1.lattice.matrix
l2 = struc_2.lattice.matrix

new_matrix = [[l2[0,0], l2[0,1], 0],
	[l2[1,0], l2[1,1], 0],
	[0, 0, 25]]

struc_1_coords = struc_1.cart_coords
struc_2_coords = struc_2.cart_coords

struc_1_species = [site.species for site in struc_1.sites]
struc_2_species = [site.species for site in struc_2.sites]

# Positions struc 1 values to be directly below z = 0
struc_1_coords[:,2] -= struc_1_coords[:,2].max()

# Positions struc 2 values to be directly at z = argv[3]
struc_2_coords[:,2] -= (struc_2_coords[:,2].min() - float(argv[3]))

hetero_coords = np.append(struc_1_coords, struc_2_coords, axis=0)

new_struc = Structure(new_matrix, species=struc_1_species + struc_2_species,
	coords= hetero_coords, coords_are_cartesian=True)

c_vec = new_struc.lattice.matrix[2]

struc_zs = new_struc.cart_coords[:, 2]
struc_dist = struc_zs.max() - struc_zs.min()
prev_centre = struc_zs.min() + (struc_dist / 2)

translation = c_vec / 2 - np.array([0, 0, prev_centre])

new_struc.translate_sites(np.arange(new_struc.cart_coords.shape[0]), vector=translation, frac_coords=False, to_unit_cell=True)

new_struc.get_sorted_structure().to(filename='POSCAR_hetero')



