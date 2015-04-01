
import Sire.Stream

from Sire.Mol import *

mol = Sire.Stream.load("../io/ligand.s3")

def printMatrix(matrix):
    for i in range(0,len(matrix)):
        print("%6d:" % i, end=" ")
        for j in range(0,len(matrix[i])):
            print("%d" % matrix[i][j], end=" ")

        print("\n", end="")

def assertSymmetrical(matrix):
    for i in range(0,len(matrix)):
        for j in range(0,len(matrix[i])):
            assert( matrix[i][j] == matrix[j][i] )

def test_matrix(verbose=False):

    connectivity = mol.property("connectivity")

    matrix = connectivity.getBondMatrix(0)

    if verbose:
        print("\nConnectivity 0")
        printMatrix(matrix)

    for i in range(0,len(matrix)):
        for j in range(0,len(matrix[i])):
            assert( matrix[i][j] )
 
    matrix1 = connectivity.getBondMatrix(1,1)

    if verbose:
        print("\nConnectivity 1")
        printMatrix(matrix1)

    for i in range(0,len(matrix1)):
        for j in range(0,len(matrix1[i])):
            assert( matrix1[i][j] == (i==j) )

    matrix2 = connectivity.getBondMatrix(2,2)

    if verbose:
        print("\nConnectivity 2")
        printMatrix(matrix2)

    assertSymmetrical(matrix2)

    for i in range(0,len(matrix2)):
        for j in range(0,len(matrix2[i])):
            assert( connectivity.areBonded( AtomIdx(i), AtomIdx(j) ) == matrix2[i][j] )

    matrix3 = connectivity.getBondMatrix(3,3)

    if verbose:
        print("\nConnectivity 3")
        printMatrix(matrix3)

    assertSymmetrical(matrix3)

    for i in range(0,len(matrix3)):
        for j in range(0,len(matrix3[i])):
            assert( connectivity.areAngled( AtomIdx(i), AtomIdx(j) ) == matrix3[i][j] )

    matrix4 = connectivity.getBondMatrix(4,4)

    if verbose:
        print("\nConnectivity 4")
        printMatrix(matrix4)

    assertSymmetrical(matrix4)

    for i in range(0,len(matrix4)):
        for j in range(0,len(matrix4[i])):
            assert( connectivity.areDihedraled( AtomIdx(i), AtomIdx(j) ) == matrix4[i][j] )

    matrix = connectivity.getBondMatrix(4)

    if verbose:
        print("\nConnectivity 1-4")
        printMatrix(matrix)

    assertSymmetrical(matrix)

    for i in range(0,len(matrix)):
        for j in range(0,len(matrix)):
            assert( matrix[i][j] == (matrix1[i][j] or matrix2[i][j] or matrix3[i][j] or matrix4[i][j]) )

if (__name__ == "__main__"):
    test_matrix(True)

