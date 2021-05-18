#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

struct matrix {
	int ncols;
	int nrows;
	double* mat;
};

void readMatrix(struct matrix* m, FILE* file);
void storeMatrix(double *squared_matrix, int n, FILE* file);
double determinant(double *m, int n);
void getCofactor(double *m, double *cofact, int p, int q, int n);
void adjoint(struct matrix *m, double *adj);
void inverse(struct matrix *m, double *inv, double det);

/*
	Reads a matrix from a file and stores it into the appropriate structure.
*/
void readMatrix(struct matrix* m, FILE* file) {
	int i, j;

	m->mat = (double*)malloc(m->ncols * m->nrows * sizeof(double));

	for (i = 0; i < m->nrows; i++) {
		for (j = 0; j < m->ncols; j++) {
			fscanf(file, "%lf", &m->mat[i * m->ncols + j]);
		}
	}
}

/*
	Stores a matrix into the file passed as argument
*/
void storeMatrix(double *squared_matrix, int n, FILE* file) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fprintf(file, "%lf ", squared_matrix[i * n + j]);
		}
		fprintf(file, "\n");
	}
}

double determinant(double *m, int n){
	double d = 0; // Initialize result
	
	//  Base case : if matrix contains single element
	if (n == 1){
		return m[0];
	}

	int sign = 1, f;  // To store sign multiplier
	double *cofact = (double*)malloc(n * n * sizeof(double));  // cofact is used to store cofactors of m

	// Iterate for each element of first row
	for (f = 0; f < n; f++){
		// Getting Cofactor of A[0][f]
		getCofactor(m, cofact, 0, f, n);

		d += sign * m[f] * determinant(cofact, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}
	
	free(cofact);
	return d;
}

void getCofactor(double *m, double *cofact, int p, int q, int n){
	int i = 0, j = 0, row, col;

	// Looping for each element of the matrix
	for (row = 0; row < n; row++){
		for (col = 0; col < n; col++){
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q){
				cofact[i * n + (j++)] = m[row * n + col];
				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1){
					j = 0;
					i++;
				}
			}
		}
	}
}

void adjoint(struct matrix *m, double *adj){
	if (m->nrows == 1) {
		adj[0] = 1;
		return;
	}

	int sign = 1, i, j;
	double *cofact = (double*)malloc(m->ncols * m->nrows * sizeof(double));  // cofact is used to store cofactors of m.mat

	for (i=0; i<m->nrows; i++){
		for (j=0; j<m->ncols; j++){
			// Get cofactor of A[i][j]
			getCofactor(m->mat, cofact, i, j, m->nrows);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i+j)%2==0)? 1: -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j * m->ncols + i] = (sign)*(determinant(cofact, m->nrows-1));
		}
	}
	free(cofact);
}

void inverse(struct matrix *m, double *inv, double det){
	int i, j;

	// Find adjoint
	double *adj = (double*)malloc(m->ncols * m->nrows * sizeof(double));
	adjoint(m, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (i=0; i<m->nrows; i++){
		for (j=0; j<m->ncols; j++){
			inv[i*m->ncols +j] = adj[i*m->ncols + j]/(det);
		}
	}
	free(adj);
}

int main(int argc, char* argv[]) {
	if(argc != 2) { //Checking parameters: 1.mat_inv.exe 2.matrix
		printf("Parameters error.\n");
		exit(1);
	}

	FILE *mat, *resultFile;
	clock_t t;
	struct matrix m;
	double det;

	mat = fopen(argv[1], "r");
	fscanf(mat, "%d %d", &m.nrows, &m.ncols);
	readMatrix(&m, mat);

	det = determinant(m.mat, m.nrows);
	printf("Determinant of the matrix: %f \n", det);

	/* Checking if it is possible to perform the matrix inversion */
	if (det == 0 || (m.nrows != m.ncols)) {
		printf("ERROR: It is not possible to compute the inversion: determinant is equal to 0 or the matrix is not squared\n");
		fclose(mat);
		fclose(resultFile);
		free(m.mat);
		exit(1);
	}

	double *inv =(double*)malloc(m.ncols * m.nrows * sizeof(double));

	t = clock();
	inverse(&m, inv, det);
	t = clock() - t;

	resultFile = fopen("inverse.txt", "w");
	storeMatrix(inv, m.nrows, resultFile);

	printf("\nElapsed time: %lf seconds\n", ((double)t) / CLOCKS_PER_SEC);

	fclose(mat);
	fclose(resultFile);
	free(inv);
	free(m.mat);

	return 0;
}
