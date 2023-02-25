import numpy as np
import math
import matplotlib.pyplot as plt


def lattice_transformations(a1, b1, a2, b2):

	"""

	Calculates the transformations and re-orientations required to transform from a1, b1 to a2, b2

	a(i) and b(i) are vectors that form 2d lattice of a single crystal
	| ai[1] ai[2] |
	| bi[1] bi[2] |

	Intended formula:
	  L2    =    M      *   L1 - is a square matrix so might be able to take inverse
	| a2 | = | i j | * | a1 |
	| b2 |   | 0 m |   | b1 |

	Where: 
	i * m = n
	i, m > 0
	0 <= j <= m - 1

	The n value above comes from:
	n = r1/r2 = area1/area2
	r1 = int, number of cells in supercell lattice1
	r2 = int, number of cells in supercell lattice2

	Output:
	M =
	| x1 x2 |
	| 0  x3 |

	"""


	# Define the initial and final lattices
	lattice1 = np.array[a1, b1]
	lattice2 = np.array[a2, b2]

	# Primative lattice transformations
	transformations = np.zeros((2,2))

	# TODO: calculate areas of each lattice and find numerical ratios (EQTN 3.1)

	# Lattice 1 area calculations
	area1 = np.linalg.norm(np.cross(a1, b1))

	# Lattice 2 area calculations
	area2 = np.linalg.norm(np.cross(a2, b2))

	# n calculation
	ratio = area2 / area1

	# QUESTION: is there a better way to store the nvals
	nVals = []

	for i in range(1,50):
		for j in range(1,50):
			if (abs(i/j - ratio)/ratio <= 0.01):
				# If ints i and j create a ratio with insignificant mismatch, append to lists
				nVals.append(i)
				nVals.append(j)

	print("N vals: ")
	print(nVals)

	x1Vals = []
	x2Vals = []
	x3Vals = []

	# Calculations: Check if the current M transforms L1 to L2
	for n in nVals:
		for x1 in range(0,n1):
			for x2 in range(0,n1):
				for x3 in ramge(0,n1):
					if (x1 * x3 == n1) and (j <= m - 1):
						m[0,0] = x1
						m[1,0] = x2
						m[1,1] = x3
						if np.dot(transformations,lattice1) == lattice2:
							x1Vals.append(x1)
							x2Vals.append(x2)
							x3Vals.append(x3)

	# TODO: Figure out which set of the x1 x2 and x3 is the best transformation







								




