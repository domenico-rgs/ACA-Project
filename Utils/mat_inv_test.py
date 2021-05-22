import sys
import numpy

def read_matrix(filename):
    mat = []
    with open(filename) as f:
        mat = [line.split() for line in f]
    # Converting strings to integers 
    for el in mat:
            for i in range(0, len(el)):
                el[i] = float(el[i])
    return mat

def store_matrix(matrix):
    result = open('inverse.txt', 'w')
    result.write(str(len(matrix)) + " " + str(len(matrix[0])) + "\n")
    for el in matrix:
        for number in el:
            result.write(str(number) + " ")
        result.write("\n")

# Checking arguments
if len(sys.argv) != 2:
    print("No argument provided, please provide a matrix.\n")
    sys.exit(1);

# Reading the two matrix
mat = read_matrix(sys.argv[1])

# Extracting first row of the matrix, containing number of rows and columns
rows_cols = mat.pop(0)

# Computing the determinant of the matrix
det = numpy.linalg.det(mat)

# Checking if it is possible to perform inversion
if (rows_cols[0] != rows_cols[1]):
    print("It is not possible to compute matrix inversion. Check matrix number of rows and cols.")
    sys.exit(1);
if (det == 0):
    print("It is not possible to compute matrix inversion. Determinant is equal to 0.")
    sys.exit(1);

# Performing inversion
inverse = numpy.linalg.inv(mat)

# Storing result matrix into a file
store_matrix(inverse)
print("Inverse matrix has been stored into inverse.txt")
