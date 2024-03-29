import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sys import argv


def calculateAreaRatio(lattice1, lattice2):
    """
    STEP 1
    Calculates the ratio between the areas of the 2 selected crystals
    Equation 2.1 of Lattice Match: An Application to heteroepitaxy
    Example shown in equation 3.1

    Lattice 1 is smaller than Lattice 2
    """

    # Lattice 1 area calculations
    area1 = np.linalg.norm(abs(np.cross(lattice1[0], lattice1[1])))

    # Lattice 2 area calculations
    area2 = np.linalg.norm(abs(np.cross(lattice2[0], lattice2[1])))

    # n calculation
    ratio = abs(area2 / area1)

    return ratio


def returnRatio(ratio, N, maxErr):
    if (ratio > 1):
        a, b = rationalizeRatio(1 / ratio, N, maxErr)
        return b, a
    elif(ratio == 1):
        return 1, 1
    else:
        return rationalizeRatio(ratio, N, maxErr)


def rationalizeRatio(ratio, N, maxErr):
    """
    STEP 2
    Calculates a set of numbers that form a rational number with a ratio
    Equation 2.2 of Lattice Match: An Application to heteroepitaxy
    Example shown in equation 3.1
    Taken from www.johndcook.com/blog/2010/10/20/best-rational-approximation/

    Intakes number 0 < x < 1
    Maximum denominator = N

    TODO: Improve rationalizeRatio algorithm
    Have rationalizeRatio test ratios within +-5%
    """

    assert ratio < 1

    a = 0
    b = 1
    c = 1
    d = 1

    while (b <= N and d <= N):
        mediant = float(a + c) / (b + d)
        # checks if it is within the user input accepted error
        if (abs(ratio - mediant) / ratio <= maxErr):
            if b + d <= N:
                return a + c, b + d
            elif d > b:
                return c, d
            else:
                return a, b
        elif ratio > mediant:
            a, b = a + c, b + d
        else:
            c, d = a + c, b + d

    if (b > N):
        return c, d
    else:
        return a, b


def calculateIndividualMVals(n):

    """
    STEP 3
    Determines the potential M matrices that represents the transformations
    to get from L1 to L2
    Equation 2.3 of Lattice Match: An Application to heteroepitaxy
    Example shown in 3.2 and 3.3

    Formula:
    L2 = M * L1
    | a2 | = | i j | * | a1 |
    | b2 | = | 0 m | * | b1 |

    Where:
    i * m = n
    i, m > 0
    0 <= j <= m - 1
    n = integer number obtained from returnRatio()

    x1 = i
    x2 = j
    x3 = m
    """

    solutions = []

    for x3 in range(0, n + 1):
        for x1 in range(0, n + 1):
            for x2 in range(0, n + 1):
                if ((x1 * x3 == n) and (x2 <= (x3 - 1))):
                    solutions.append([[x1, x2], [0, x3]])

    return solutions


def determineBestMatrix(lattice1, lattice2, mMatrices, nMatrices):

    """
    STEP 4
    Tests sets of testMatrices to determine which matrixes result
    in the smallest difference between lattice1 and lattice2
    Determines which angle has the smallest difference in superlattices
    """

    # Large initial error to compare against
    smallestErr = 1000000

    # Iterates to test each pair of m and n matrices for lowest error
    for m in mMatrices:
        for n in nMatrices:
            for theta in np.linspace(0, 360, 3600):
                # TODO: Incorporate algorithm to test angles
                # Allow it to choose any angle 0-360, figure out the optimal solution

                # Represents the crystal after
                # lattices have undergone transformations
                transformedL1 = np.dot(lattice1, m)
                transformedL2 = np.dot(lattice2, n)

                radianTheta = np.radians(theta)

                # Temp matrix to test counterclockwise rotation
                rotation = [[np.cos(radianTheta), -1 * np.sin(radianTheta)],
                    [np.sin(radianTheta), np.cos(radianTheta)]]

                # Applies rotation matrix to L1
                transformedL1 = np.dot(transformedL1, rotation)

                diffMatrix = (transformedL2 - transformedL1) ** 2

                diff1 = diffMatrix[0][0]
                diff2 = diffMatrix[0][1]
                diff3 = diffMatrix[1][0]
                diff4 = diffMatrix[1][1]

                # Calculation of the RMS between the two matrices
                diff = np.sqrt((diff1 + diff2 + diff3 + diff4)/4)

                # Updates the best set of transformations
                # if current set has smaller diff
                if (diff < smallestErr):
                    smallestErr = diff
                    bestM = m
                    bestN = n
                    bestR = theta

    return [bestM, bestN, bestR, smallestErr]


def graphSuperLattices(lattice1, lattice2, acceptableMatrices, nVals):

    """
    STEP 5
    Graphs the superlattices imposed on each other
    """

    m = acceptableMatrices[0]
    n = acceptableMatrices[1]
    theta = acceptableMatrices[2]

    # Calculates the rotation matrix
    radianTheta = np.radians(theta)
    rotation = [[np.cos(radianTheta), -1 * np.sin(radianTheta)],
                [np.sin(radianTheta), np.cos(radianTheta)]]

    # Calculates the superlattices that will be plotted
    transformedL1 = np.dot(np.dot(lattice1, m), rotation)
    transformedL2 = np.dot(lattice2, n)

    # Sets the points of the parallelogram created by transformedL1
    p12 = [transformedL1[0][0], transformedL1[0][1]]
    p13 = [transformedL1[0][0] + transformedL1[1][0],
            transformedL1[0][1] + transformedL1[1][1]]
    p14 = [transformedL1[1][0], transformedL1[1][1]]

    # Creates lists of the x and y points of L1
    l1XPoints = [0, p12[0], p13[0], p14[0]]
    l1YPoints = [0, p12[1], p13[1], p14[1]]

    # Sets the points of the parallelogram created by transformedL2
    p22 = [transformedL2[0][0], transformedL2[0][1]]
    p23 = [transformedL2[0][0] + transformedL2[1][0],
            transformedL2[0][1] + transformedL2[1][1]]
    p24 = [transformedL2[1][0], transformedL2[1][1]]

    # Creates lists of the x and y points of L2
    l2XPoints = [0, p22[0], p23[0], p24[0]]
    l2YPoints = [0, p22[1], p23[1], p24[1]]

    # Sets the minimum and maximum x boundaries for graph
    allXPoints = l1XPoints + l2XPoints
    minX = min(allXPoints) * 1.1
    maxX = max(allXPoints) * 1.1

    # Sets the minimum and maximum y boundaries for graph
    allYPoints = l1YPoints + l2YPoints
    minY = min(allYPoints) * 1.1
    maxY = max(allYPoints) * 1.1

    # If the min x or Y val is 0, add margin to graph
    if (minX == 0):
        minX = minX - 0.1 * maxX

    if (minY == 0):
        minY = minY - 0.1 * maxY

    # Set up superlattice graph
    fig, ax = plt.subplots()
    plt.grid(True)
    plt.title("Superlattice Mismatch")
    ax.set_xlabel('X axis (Angstroms)')
    ax.set_ylabel('Y axis (Angstroms)')

    # Set graph boundaries
    plt.xlim([minX, maxX])
    plt.ylim([minY, maxY])

    # Graph superlattice 1
    ax.add_patch(patches.Polygon(xy=list(zip(l1XPoints, l1YPoints)),
        fill=False, color='red', linewidth=2))

    # Graph superlattice 2
    ax.add_patch(patches.Polygon(xy=list(zip(l2XPoints, l2YPoints)),
        fill=False, color='blue', linewidth=2))

    plt.show()


def latticeTransformations(lattice1, lattice2, maxN, maxErr):

    """
    Calculates the transformations and re-orientations
    required to transform from a1, b1 to a2, b2

    a(i) and b(i) are vectors that form 2d lattice of a single crystal
    | ai[1] ai[2] |
    | bi[1] bi[2] |

    Intended formula:
    L2 = M * L1 - is a square matrix so might be able to take inverse
    | a2 | = | i j | * | a1 |
    | b2 | = | 0 m | * | b1 |

    Where:
    i * m = n
    i, m > 0
    0 <= j <= m - 1

    The n value above comes from:
    n = r1/r2 = area1/area2
    r1 = int, number of cells in supercell lattice1
    r2 = int, number of cells in supercell lattice2

    Output:
    M = [[x1, x2], [0, x3]]
    """

    # List of sets of n values that create the area ratios
    nVals = []

    # List of values [[x1, 0], [x2, x3]] that compose the M transformation matrix
    mMatrices = []

    acceptableMatrices = []

    # STEP 1
    # Calculates the area ratio between lattice 1 and lattice 2
    ratio = calculateAreaRatio(lattice1, lattice2)

    # STEP 2
    # Calculates the integer ratio of areas
    # TO FIX: 
    nVals = returnRatio(ratio, maxN, maxErr)

    # STEP 3
    # Calculate the possible M matrices
    # mMatrices: matrices to multiply lattice2 by
    # nMatrices: matrices to multiply lattice1 by
    mMatrices = calculateIndividualMVals(nVals[0])
    nMatrices = calculateIndividualMVals(nVals[1])

    # STEP 4
    # Determine which matrix set has the lowest error
    acceptableMatrices = determineBestMatrix(lattice1, 
        lattice2, mMatrices, nMatrices)

    # STEP 5
    # Graph the acceptableMatrices results
    graphSuperLattices(lattice1, lattice2, acceptableMatrices, nVals)

    return acceptableMatrices

