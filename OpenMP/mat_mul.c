#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void readMatrix(double** m, FILE* file, int rows, int cols);
void printMatrix(double** m, FILE* file, int rows, int cols);
void matrixMul(double** m1, double** m2, double** m3, int m, int p);

/*
 * The file which contains a matrix has in its first row the dimensions
 * then using fscanf each element of the matrix is stored on the memory allocated dynamically
*/
void readMatrix(double** m, FILE* file, int rows, int cols) {
	int i, j;

	for(i = 0; i < rows; i++) {
		m[i] = (double*)malloc(cols*sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			fscanf(file, "%lf", &m[i][j]);
		}
	}
}

/* The opposite operation of readMatrix. Stores a matrix into a file, element by element */
void printMatrix(double** m, FILE* file, int rows, int cols) {
	int i, j;

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			fprintf(file, "%lf ", m[i][j]);
		}
		fprintf(file, "\n");
	}
}

/*
 * Performs the multiplication operation between the matrices m1 and m2.
 * The result will be stored in the matrix m3.
 * The algorithm is practically the one that can be found here: https://en.wikipedia.org/wiki/Matrix_multiplication#Definition
*/
void matrixMul(double** m1, double** m2, double** m3, int m, int p) {
	int i, j, k;
	
	for(i = 0; i < m; i++) {
		m3[i] = (double*)malloc(p * sizeof(double));
		memset(m3[i], 0,  p * sizeof(double)); 
	}
	
	#pragma omp parallel shared(m1, m2, m3) private(i, j, k)
	{
		#pragma omp for schedule (dynamic) //the iterations of the first for loop are shared among threads, this means that each thread calculates a row of the resulting matrix
		for (i = 0; i < m; i++) {
			for (j = 0; j < p; j++) {
				for (k = 0; k < p; k++) { 
					m3[i][j] += m1[i][k] * m2[k][j]; //there are no data races since each thread access independetly to only one m3 element at a time
				}	
			}
		}
	}
}

int main(int argc, char* argv[]) {
	if(argc != 3){ //1- exe name, 2- mat1.txt, 3- mat2.txt
		printf("Parameter error.\n");
		exit(1);
	}

	FILE *mat1, *mat2, *resultFile;
	double t;
	int m, n1, n2, p, i;

	mat1 = fopen(argv[1], "r");
	mat2 = fopen(argv[2], "r");
	fscanf(mat1, "%d %d", &m, &n1);
	fscanf(mat2, "%d %d", &n2, &p);

	/* Multiplication is permitted if m1 is m x n and m2 is n x p, m1 must have the same number of column of the rows of m2 matrix */
	if(n1 != n2) {
		printf("It is not possible to do matrix multiplication. Check matrix number of rows and cols.\n");
		fclose(mat1);
		fclose(mat2);
		exit(1);
	}
	
	double **m1 = (double **)malloc(m*sizeof(double*));
	double **m2 = (double **)malloc(n2*sizeof(double*));
	double **m3 = (double **)malloc(m*sizeof(double*));
	
	readMatrix(m1, mat1, m, n1);
	readMatrix(m2, mat2, n2, p);
	
	t = omp_get_wtime();
	matrixMul(m1, m2, m3, m, p);
	t = omp_get_wtime() - t; //total time spent in matrixMul (wall clock time)

	resultFile = fopen("result.txt", "w");
	printMatrix(m3, resultFile, m, p);

	printf("Elapsed time: %.5f seconds\n", t);

	fclose(mat1);
	fclose(mat2);
	fclose(resultFile);

	for(i = 0; i < m; i++) {
		free(m1[i]);
		free(m3[i]);
	}
	for(i = 0; i < n2; i++) {
		free(m2[i]);
	}
	free(m1);
	free(m2);
	free(m3);
	return 0;
}
