import numpy as np
import math
#import matplotlib.pyplot as plt


def calculateAreaRatio(lattice1, lattice2):
	""" 
	STEP 1
	Calculates the ratio between the areas of the 2 selected crystals
	Equation 2.1 of Lattice Match: An Application to heteroepitaxy
	Example shown in equation 3.1

	Lattice 1 is smaller than Lattice 2
	"""

	# Lattice 1 area calculations
	#np.linalg.det() - google this to figure out how to get area and take abs of it
	area1 = np.linalg.norm(np.cross(lattice1[0], lattice1[1]))

	# Lattice 2 area calculations 
	area2 = np.linalg.norm(np.cross(lattice2[0], lattice2[1]))

	# n calculation
	ratio = area2 / area1

	return ratio

# -------------------------------------------------------------------------------

def rationalizeRatio(ratio, N):
	"""
	STEP 2
	Calculates a set of numbers that form a rational number with a ratio
	Equation 2.2 of Lattice Match: An Application to heteroepitaxy
	Example shown in equation 3.1
	Taken from # https://www.johndcook.com/blog/2010/10/20/best-rational-approximation/

	Intakes number 0 < x < 1
	Maximum denominator = N
	"""

	a = 0
	b = 1
	c = 1
	d = 1

	# Flip ratio if it is greater than 1 to fit rationalizeRatio
	if (ratio > 1):
		ratio = 1 / ratio
	elif (ratio == 1):
		return 1,1

	while (b <= N and d <= N):
		mediant = float(a+c)/(b+d)
		# checks if it is within the 1% accepted error
		if (abs(ratio - mediant)/ratio <= 0.01):
			if b + d <= N:
				return a+c, b+d
			elif d > b:
				return c, d
			else:
				return a, b
		elif ratio > mediant:
			a, b = a+c, b+d
		else:
			c, d = a+c, b+d

	if (b > N):
		return c, d
	else:
		return a, b

# -------------------------------------------------------------------------------

def calculateIndividualMVals(n):

	"""
	STEP 3
	Determines the potential M matrices that represents the transformations 
	to get from L1 to L2
	Equation 2.3 of Lattice Match: An Application to heteroepitaxy
	Example shown in  3.2 and 3.3

	Formula: 
	  L2    =    M      *   L1 - is a square matrix so might be able to take inverse
	| a2 | = | i j | * | a1 |
	| b2 |   | 0 m |   | b1 |

	Where: 
	i * m = n
	i, m > 0
	0 <= j <= m - 1
	n = integer number obtained from rationalizeRatio()

	x1 = i
	x2 = j
	x3 = m

	"""

	m = [[0, 0], [0, 0]]

	solutions = []

	for x3 in range(0, n + 1):
		for x1 in range(0, n + 1):
			for x2 in range(0, n + 1):
				if ((x1 * x3 == n) and (x2 <= (x3 - 1))):
					solutions.append([[x1, x2], [0,x3]])

					# print("x1 = " + str(x1))
					# print("x2 = " + str(x2))
					# print("x3 = " + str(x3))

	return solutions

# -------------------------------------------------------------------------------

def determineBestMatrix(lattice1, lattice2, mMatrices, nMatrices):

	"""
	STEP 4
	Tests individual testMatrices to see if they are within the acceptable 1% error 
	"""

	# Variables to store the best M, N, and error set
	bestM = []
	bestN = []
	smallestErr = 100000000

	# Iterates to test each pair of m and n matrices for lowest error
	for m in mMatrices:
		for n in nMatrices:
			# Represents the crystal after lattices have undergone transformations to be close 
			transformedL1 = np.dot(lattice1, n)
			transformedL2 = np.dot(lattice2, m)

			# Calculates the root mean square difference between matrices
			diff = ((transformedL2[0][0] - transformedL1[0][0]) ** 2) + \
			((transformedL2[0][1] - transformedL1[0][1]) ** 2) + \
			((transformedL2[1][0] - transformedL1[1][0]) ** 2) + \
			((transformedL2[1][1] - transformedL1[1][1]) ** 2)

			# Updates the best set of transformations if current set has smaller diff
			if diff < smallestErr:
				bestM = m
				bestN = n
				smallestErr = diff

	return [bestM, bestN, smallestErr]

# -------------------------------------------------------------------------------

def lattice_transformations(lattice1, lattice2):

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

	# Primative lattice transformations
	transformations = np.zeros((2,2))

	# List of sets of n values that create the area ratios
	nVals = []

	# List of values [[x1, 0], [x2, x3]] that compose the M transformation matrix
	mMatrices = []

	acceptableMatrices = []

	# STEP 1
	# Calculates the area ratio between lattice 1 and lattice 2
	ratio = calculateAreaRatio(lattice1, lattice2)

	# STEP 2
	# Calculates the integer ratio with 
	nVals = rationalizeRatio(ratio, 1000)


	# STEP 3
	# Calculate the possible M matrices given 
	# mMatrices: matrices to multiply lattice2 by
	# nMatrices: matrices to multiply lattice1 by
	mMatrices = calculateIndividualMVals(nVals[0])
	nMatrices = calculateIndividualMVals(nVals[1])

	# STEP 4
	# Determine which matrix set has the lowest error
	acceptableMatrices = determineBestMatrix(lattice1, lattice2, mMatrices, nMatrices)

	return acceptableMatrices
	
	# TODO: Figure out which set of the x1 x2 and x3 is the best transformation

# Future final test case for overall script

# -------------------------------------------------------------------------------

# # Unit tests for calculateAreaRatio()
# print("Testing calculateAreaRatio()")
# print()

# Test of 1 to 1
# print("Identical Lattices: " + str(calculateAreaRatio([[1,0], [0,1]], [[1,0], [0,1]])))
# print("Flipped Lattices: " + str(calculateAreaRatio([[1,0], [0,1]], [[0,1], [1,0]])))


# # Test 9/4 = 5.6025
# print("Testing 9/4: " + str(calculateAreaRatio([[4,0], [0,4]], [[9,0], [0,9]])))

# # Testing non-integer ratios
# print("Testing 9.1/4: " + str(calculateAreaRatio([[4,0], [0,4]], [[9.1,0], [0,9.1]])))

# # GaAs: a = 5.653
# # CdTe: a = 6.481
# # Ratio should be 1.314
# print("GaAs and CdTe ratio: " + str(calculateAreaRatio([[5.653,0], [0,5.653]], [[6.481,0], [0,6.481]])))

# print("----------------")

# Unit tests for rationalizeRatio()
# print("Testing rationalizeRatio()")
# print()

# Testing known number ratios
# print("1/2: " + str(rationalizeRatio(0.500, 1000)))
# print("1/3: " + str(rationalizeRatio(0.3333333, 1000)))
# print("2/3: " + str(rationalizeRatio(0.6666667, 1000)))
# print("1/4: " + str(rationalizeRatio(0.25, 1000)))
# print("1/1: " + str(rationalizeRatio(1.000, 1000)))

# print()

# # Testing GaAs and CdTe
# print("Testing GaAs and CdTe: " + str(rationalizeRatio(calculateAreaRatio([[5.653,0], [0,5.653]], [[6.481,0], [0,6.481]]), 1000)))
# print("Testing 1.314: " + str(rationalizeRatio(1.314395525, 1000)))
# print("Testing 0.760: " + str(rationalizeRatio(0.7608059984, 1000)))

# print()

# # Testing known irrational numbers
# print("Testing pi: " + str(rationalizeRatio(math.pi , 1000)))
# print("Testing e:" + str(rationalizeRatio(math.e, 1000)))
# print("Testing sqrt(2):" + str(rationalizeRatio(math.sqrt(2), 1000)))

# print("----------------")

# Unit tests for calculateM()
# print("Testing calculateIndividualMVals()")
# print()

# print("n = 0: " + str(calculateIndividualMVals(0)))
# print()

# print("n = 1: " + str(calculateIndividualMVals(1)))
# print()

# print("n = 2: " + str(calculateIndividualMVals(2)))
# print()

# print("n = 3: " + str(calculateIndividualMVals(3)))
# print()

# print("n = 4: " + str(calculateIndividualMVals(4)))
# print()

# print("n = 5: " + str(calculateIndividualMVals(5)))
# print()

# print("n = 6: " + str(calculateIndividualMVals(6)))
# print()

# print("----------------")

# print("Testing calculateAllMVals()")
# print()

# print("Set 0,1: " + str(calculateAllMVals([0,1])))

# print("Set 1,2: " + str(calculateAllMVals([1,2])))

# # The n values used for CdTe GaAs example, provides same outpouts
# print("Set 3,4: " + str(calculateAllMVals([3,4])))

# print("----------------")
# print()

# print("Testing calculatePercentError()")
# print()

# Standard test with simple numbers easy to calculate with
# print(calculatePercentError([[1,0], [0,1]], [[2,0], [0,2]], [[1, 0], [0, 1]], [[2, 0], [0, 2]]))

# Testing lattice matrix 
# print(calculatePercentError([[1,0], [0,1]], [[0,1], [1,0]], [[1, 0], [0, 1]]))

# Testing a case with an acceptable erro
# print(calculatePercentError([[1,0], [0,1]], [[1.001,0], [0,1.001]], [[1, 0], [0, 1]]))

# print("----------------")
# print()

print("Testing lattice_transformations()")
print()

print("Testing identical lattices: " + str(lattice_transformations([[1,0], [0,1]], [[1,0], [0,1]])))

print("Testing 2x Lattice: " + str(lattice_transformations([[1,0], [0,1]], [[2,0], [0,2]])))

# print("Testing 2x Lattice with rotation " + str(lattice_transformations([[1,0], [0,1]], [[0,2], [2,0]])))

print("Testing 3x Lattices: " + str(lattice_transformations([[1,0], [0,1]], [[3,0], [0,3]])))

# print("Testing 3x Lattice with rotation: " + str(lattice_transformations([[1,0], [0,1]], [[0,3], [3,0]])))

# print("Testing 4x Lattices: " + str(lattice_transformations([[1,0], [0,1]], [[4,0], [0,4]])))

# print("Testing Non-rectangular Lattice2: " + str(lattice_transformations([[1,0], [0,4]], [[1,0], [0,4]])))

# print("Testing 4x Lattice with rotation: " + str(lattice_transformations([[1,0], [0,1]], [[0,4], [4,0]])))

# print("Testing CdTe and GaAs: " + str(lattice_transformations([[5.653,0], [0,5.653]], [[6.481,0], [0,6.481]])))

# http://www.2dmatpedia.org/2dmaterials/doc/2dm-2997 HgBrN
# print("Testing 1 Angstrom to HgBrN: " + str(lattice_transformations([[1,0], [0,1]], [[4.02,0], [0,4.46]])))

# http://www.2dmatpedia.org/2dmaterials/doc/2dm-2994 Ga: 2.66
# http://www.2dmatpedia.org/2dmaterials/doc/2dm-2998 LiMg: 3.18
# http://www.2dmatpedia.org/2dmaterials/doc/2dm-3000 VTe2: 3.65
# http://www.2dmatpedia.org/2dmaterials/doc/2dm-2995 Sb2Te3: 4.32

print("Testing Ga to LiMg: " + str(lattice_transformations([[2.66,0], [0,2.66]], [[3.18,0], [0,3.18]])))

print("Testing Ga to Sb2Te3: " + str(lattice_transformations([[2.66,0], [0,2.66]], [[4.32,0], [0,4.32]])))

print("Testing LiMg to Sb2Te3: " + str(lattice_transformations([[3.18,0], [0,3.18]], [[4.32,0], [0,4.32]])))









