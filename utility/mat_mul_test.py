import sys

def read_matrix(filename):
    mat = []
    with open(filename) as f:
        mat = [line.split() for line in f]
    # Converting strings to integers 
    for el in mat:
            for i in range(0, len(el)):
                el[i] = float(el[i])
    return mat

def store_matrix(mat):
    result = open('result.txt', 'w')
    result.write(str(len(mat)) + " " + str(len(mat[0])) + "\n")
    for el in mat:
        for number in el:
            result.write(str(number) + " ")
        result.write("\n")

# Checking arguments
if len(sys.argv) != 3:
    print("No argument provided, please provide two matrix files.\n")
    sys.exit(1);

# Reading the two matrix
mat1 = read_matrix(sys.argv[1])
mat2 = read_matrix(sys.argv[2])

# Extracting first row of the two matrix, containing number of rows and columns
rows_cols1 = mat1.pop(0)
rows_cols2 = mat2.pop(0)

# Checking if it is possible to perform multiplication
if (rows_cols1[1] != rows_cols2[0]):
    print("It is not possible to do matrix multiplication. Check matrix number of rows and cols.")
    sys.exit(1);

# Performing multiplication
mat_result = [[sum(m1 * m2 for m1, m2 in zip(mat1_row, mat2_col)) 
                        for mat2_col in zip(*mat2)]
                                for mat1_row in mat1]

# Storing result matrix into a file
store_matrix(mat_result)
print("The result matrix has been stored into result.txt")
