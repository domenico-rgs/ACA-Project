/* https://en.wikipedia.org/wiki/LU_decomposition#C_code_example */

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
double LUPDeterminant(double *m, int *P, int n);
void LUPInvert(double *m, int *P, int n, double *inverseM);
void LUPSolve(double *m, int *P, double *b, int n, double *x);
int LUPDecompose(double *m, int n, double Tol, int *P);

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

/* INPUT: A,P filled in LUPDecompose; N - dimension.
* OUTPUT: Function returns the determinant of the initial matrix
*/
double LUPDeterminant(double *m, int *P, int n){
	int i;
	double det = m[0];

	for (i=1; i<n; i++){
		det *= m[i*n+i];
	}

	return (P[n] - n) % 2 == 0 ? det : -det;
}

/* INPUT: A,P filled in LUPDecompose; N - dimension
* OUTPUT: IA is the inverse of the initial matrix
*/
void LUPInvert(double *m, int *P, int n, double *inverseM) {
	int j, i, k;

	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {
			inverseM[i*n+j] = P[i] == j ? 1.0 : 0.0;

			for (k = 0; k < i; k++){
				inverseM[i*n+j] -= m[i*n+k] * inverseM[k*n+j];
			}
		}

		for (i = n - 1; i >= 0; i--) {
			for (k = i + 1; k < n; k++){
				inverseM[i*n+j] -= m[i*n+k] * inverseM[k*n+j];
			}
			inverseM[i*n+j] /= m[i*n+i];
		}
	}
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
* OUTPUT: x - solution vector of A*x=b
*/
void LUPSolve(double *m, int *P, double *b, int n, double *x) {
	int i, k;
	for (i = 0; i < n; i++) {
		x[i] = b[P[i]];

		for (k = 0; k < i; k++){
			x[i] -= m[i*n+k] * x[k];
		}
	}

	for (i = n - 1; i >= 0; i--) {
		for (k = i + 1; k < n; k++){
			x[i] -= m[i*n+k] * x[k];
		}

		x[i] /= m[i*n+i];
	}
}

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
*        Tol - small tolerance number to detect failure when the matrix is near degenerate
* OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
*        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
*        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
*        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
*/
int LUPDecompose(double *m, int n, double Tol, int *P) {

	int i, j, k, imax;
	double maxA, absA, temp;

	for (i = 0; i <= n; i++){
		P[i] = i; //Unit permutation matrix, P[N] initialized with N
	}

	for (i = 0; i < n; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < n; k++){
			if ((absA = fabs(m[k*n+i])) > maxA) {
				maxA = absA;
				imax = k;
			}
		}

		if (maxA < Tol) return 0; //failure, matrix is degenerate

		if (imax != i) {
			//pivoting P
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			//pivoting rows of A
			for(j = 0; j < n; j++){
				temp = m[i*n+j];
				m[i*n + j] = m[imax*n+j];
				m[imax*n+j]=temp;
			}

			//counting pivots starting from n (for determinant)
			P[n]++;
		}

		for (j = i + 1; j < n; j++) {
			m[j*n+i] /= m[i*n+i];

			for (k = i + 1; k < n; k++){
				m[j*n+k] -= m[j*n+i] * m[i*n+k];
			}
		}
	}

	return 1;  //decomposition done
}

int main(int argc, char* argv[]) {
	if(argc != 2) { //Checking parameters: 1.mat_inv.exe 2.matrix
		printf("Parameters error.\n");
		exit(1);
	}

	FILE *mat, *resultFile;
	clock_t t;
	struct matrix m;

	mat = fopen(argv[1], "r");
	fscanf(mat, "%d %d", &m.nrows, &m.ncols);
	readMatrix(&m, mat);

	if (m.nrows != m.ncols) {
		printf("ERROR: It is not possible to compute the inversion: the matrix is not squared\n");
		fclose(mat);
		free(m.mat);
		exit(1);
	}

	int *P = 	(int*)malloc(m.nrows * sizeof(int) + 1 * sizeof(int));
	double *inverseM = (double*)malloc(m.ncols * m.nrows * sizeof(double));
	double *x = (double*)malloc(m.nrows * sizeof(double));
	double *b = (double*)malloc(m.nrows * sizeof(double));

	t = clock();

	double check = LUPDecompose(m.mat, m.nrows, 1E-3, P);
	
	double det = LUPDeterminant(m.mat, P, m.nrows);
	printf("\nDeterminant: %lf\n", det);
	if (det == 0.0 || check == 0) {
		printf("ERROR: It is not possible to compute the inversion\n");
		fclose(mat);
		free(P);
		free(inverseM);		
		free(x);
		free(b);
		free(m.mat);
		exit(1);
	}
	
	LUPSolve(m.mat, P, b, m.nrows, x);
	LUPInvert(m.mat, P, m.nrows, inverseM);

	t = clock() - t;

	resultFile = fopen("inverse.txt", "w");
	storeMatrix(inverseM, m.nrows, resultFile);

	printf("\nElapsed time: %lf seconds\n", ((double)t) / CLOCKS_PER_SEC);

	fclose(mat);
	fclose(resultFile);

	free(P);
	free(inverseM);
	free(x);
	free(b);
	free(m.mat);

	return 0;
}
