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
	if(ratio > 1):
		ratio = 1 / ratio

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

# Questions for Chuin Wei
def calculateIndividualMVals(n):

	"""
	STEP 3A
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

					# m = [[x1, x2], [0, x3]]

					# temp = np.dot(m, lattice1)

					# if ((np.dot(m, lattice1) == lattice2).all()):
					# 	return x1, x2, x3

	return solutions

# -------------------------------------------------------------------------------

def calculateAllMVals(lattice1, lattice2, nVals):

	"""
	STEP 3B
	Calculates M matrices for all potential N values calculated in rationalizeRatio()
	"""

	xVals = []

	for ratioPair in nVals:
		solutions1 = calculateIndividualMVals(lattice1, lattice2, ratioPair[0])
		for sol in solutions1:
			xVals.append(sol)
		solutions2 = calculateIndividualMVals(lattice1, lattice2, ratioPair[1])
		for sol in solutions2:
			xVals.append(sol)

	return xVals


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

	# List of values [x1, x2, x3] that compose the M transformation matrix
	xVals = []

	# STEP 1
	# Calculates the area ratio between lattice 1 and lattice 2
	ratio = calculateAreaRatio(lattice1, lattice2)

	# STEP 2
	# Calculates the integer ratio with 
	for i in range(1,50):
		nVals.append(rationalizeRatio(ratio, i))

	print("N vals: " + nVals)

	# STEP 3
	# Calculate the possible M matrices given 
	
	
	# TODO: Figure out which set of the x1 x2 and x3 is the best transformation

# Future final test case for overall script
#lattice_transformations(np.array([[5.653,0], [0,5.653]], [6.481,0], [0,6.481])


# -------------------------------------------------------------------------------

# # Unit tests for calculateAreaRatio()
# print("Testing calculateAreaRatio()")
# print()

# # Test of 1 to 1
# print("Identical Lattices: " + str(calculateAreaRatio([[1,0], [0,1]], [[1,0], [0,1]])))


# # Test 9/4 = 5.6025
# print("Testing 9/4: " + str(calculateAreaRatio([[4,0], [0,4]], [[9,0], [0,9]])))

# # Testing non-integer ratios
# print("Testing 9.1/4: " + str(calculateAreaRatio([[4,0], [0,4]], [[9.1,0], [0,9.1]])))

# # GaAs: a = 5.653
# # CdTe: a = 6.481
# # Ratio should be 1.314
# print("GaAs and CdTe ratio: " + str(calculateAreaRatio([[5.653,0], [0,5.653]], [[6.481,0], [0,6.481]])))

# print("----------------")

# # Unit tests for rationalizeRatio()
# print("Testing rationalizeRatio()")
# print()

# # Testing known number ratios
# print("1/2: " + str(rationalizeRatio(0.500, 1000)))
# print("1/3: " + str(rationalizeRatio(0.3333333, 1000)))
# print("2/3: " + str(rationalizeRatio(0.6666667, 1000)))
# print("1/4: " + str(rationalizeRatio(0.25, 1000)))

# # TODO: Fix when ratio = 1, algorithm does not have a numerator == denominator
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
print("Testing calculateIndividualMVals()")
print()

# Testing Simple Calculations 
# All answers checked against matrix multiplication calculator
print("n = 1: " + str(calculateIndividualMVals(1)))
print()

print("n = 2: " + str(calculateIndividualMVals(2)))
print()

print("n = 3: " + str(calculateIndividualMVals(3)))
print()

print("n = 4: " + str(calculateIndividualMVals(4)))
print()

print("n = 5: " + str(calculateIndividualMVals(5)))
print()

print("n = 6: " + str(calculateIndividualMVals(6)))
print()




