# ACA-Project
It is the project for Advanced Computer Architecture exam where we try to optimize matrix multiplication and inversion using parallel computing.

### Files
- Serial
  - mat_mul.c: receives two files with a matrix in each one as argument and multiply them giving file with the result
  - mat_inv.c: receives a file with a matrix as argument and calculates its inverse saving it in a file
- CUDA
- mat_gen.py: generates one or random matrix (depending on the number of the arguments) with a dimension specified as argument

mat_mul.c and mat_inv.c measure the time taken to perform the calculations.

### Txt matrix file example
```
3 3
4.69507 5.57633 1.72311
5.32656 6.25589 7.42191
0.86880 0.36247 5.76551
```

The first row contains the size of the matrix, number of rows and columns respectively

## How-To-Use
### mat_mul.c

1) Open a terminal in the project directory <br>

2) Compile the source code by running: <br>
``` gcc -O3  mat_mul.c ``` <br>

3) Run the program by running: <br>
``` a.exe [mat1] [mat2] ```

    Arguments:
    * mat1 (**required**): a txt file with the first matrix.
    * mat2 (**required**): second txt file with the other matrix.
