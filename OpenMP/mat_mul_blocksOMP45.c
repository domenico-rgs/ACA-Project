#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

struct matrix {
	int ncols;
	int nrows;
	double* mat;
};

#define BLOCKSIZE 64

void readMatrix(double** m, FILE* file, int rows, int cols);
void printMatrix(double** m, FILE* file, int rows, int cols);
void matrixMul(int BS, double** m1, double** m2, double** m3, int m, int p);

/*
 * Knowing the number of rows and columns,
 * it reads a matrix from a file and stores it in the appropriate structure.
*/
void readMatrix(double** m, FILE* file, int rows, int cols) {
	int i, j;

	for(i=0; i<rows; i++){
		m[i]=(double*)malloc(cols*sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			fscanf(file, "%lf", &m[i][j]);
		}
	}
}

/*
 * The opposite operation of readMatrix. Stores a matrix into the file passed as argument
*/
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
*/
void matrixMul(int BS, double** m1, double** m2, double** m3, int m, int p) {
	int i, j, k;
	int ii, jj, kk;	
	
	for(i=0; i<m; i+=BS){
		for(j=0; j<p; j+=BS){
			for(k=0; k<p; k+=BS){
				#pragma omp task depend(in: m1[i:BS][k:BS], m2[k:BS][j:BS]) depend(inout: m3[i:BS][j:BS])
				for (ii = i; ii < i+BS; ii++) {
					for (jj = j; jj < j+BS; jj++) {
						for (kk = k; kk < k+BS; kk++) { 
							m3[ii][jj] += m1[ii][kk] * m2[kk][jj];
						}	
					}
				}
			}
		}
	}
}

int main(int argc, char* argv[]) {
	if(argc != 3){ //1- exe name, 2- mat1, 3- mat2
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

	/* Multiplication is permitted if m1 is m x n and m2 is n x p */
	if(n1 != n2) {
		printf("It is not possible to do matrix multiplication. Check matrix number of rows and cols.\n");
		fclose(mat1);
		fclose(mat2);
		exit(1);
	}
	
	double ** m1 = (double **)malloc(m*sizeof(double*));
	double ** m2 = (double **)malloc(n2*sizeof(double*));
	double ** m3 = (double **)malloc(m*sizeof(double*));
	
	readMatrix(m1, mat1, m, n1);
	readMatrix(m2, mat2, n2, p);
	
	t = omp_get_wtime();
	#pragma omp parallel for private(i)
	for(i=0; i<m; i++){
		m3[i]=(double*)malloc(p*sizeof(double));
		memset(m3[i], 0,  p * sizeof(double)); 
	}
		
	#pragma omp parallel
	#pragma omp single
	matrixMul(BLOCKSIZE, m1, m2, m3, m, p);
	t = omp_get_wtime() - t; // total time spent in matrixMul

	resultFile = fopen("result.txt", "w");
	printMatrix(m3, resultFile, m, p);

	printf("Elapsed time: %.5f seconds\n", t);

	fclose(mat1);
	fclose(mat2);
	fclose(resultFile);
	for(i=0; i<m; i++){
		free(m1[i]);
		free(m3[i]);
	}
	for(i=0; i<n2; i++){
		free(m2[i]);
	}
	free(m1);
	free(m2);
	free(m3);
	return 0;
}
