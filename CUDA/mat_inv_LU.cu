/* https://en.wikipedia.org/wiki/LU_decomposition#C_code_example */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "cuda_runtime.h"
//#include "cuda_profiler_api.h"

#define THREADS 32

struct matrix {
	int ncols;
	int nrows;
	double* mat;
};

void readMatrix(struct matrix* m, FILE* file);
void storeMatrix(double *squared_matrix, int n, FILE* file);
double determinant(double *m, int n);
double LUPDeterminant(double *m, int *P, int n);
__global__ void LUPInvert(double *d_m, int *d_P, int n, double *d_inverseM);
__global__ void LUPSolve(double *d_m, int *d_P, double *d_b, int n, double *d_x);
__global__ void LUPSolve2(double *d_m, int n, double *d_x);
__global__ void LUPDecompose(double *d_m, int n, double Tol, int *d_P);

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
__global__ void LUPInvert(double *d_m, int *d_P, int n, double *d_inverseM) {
	int i = blockIdx.y*blockDim.y+threadIdx.y;
	int j = blockIdx.x*blockDim.x+threadIdx.x;
	int k;

	if((i<n) && (j<n){
		d_inverseM[i*n+j] = d_P[i] == j ? 1.0 : 0.0;

		for (k = 0; k < i; k++){
			d_inverseM[i*n+j] -= d_m[i*n+k] * d_inverseM[k*n+j];
		}

		for (i = n - 1; i >= 0; i--) {
			for (k = i + 1; k < n; k++){
				d_inverseM[i*n+j] -= d_m[i*n+k] * d_inverseM[k*n+j];
			}
			d_inverseM[i*n+j] /= d_m[i*n+i];
		}
	}
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
* OUTPUT: x - solution vector of A*x=b
*/
__global__ void LUPSolve(double *d_m, int *d_P, double *d_b, int n, double *d_x) {
	int i = blockIdx.x * THREADS + threadIdx.x;
	int k;

	if(i<n){
		d_x[i] = d_b[d_P[i]];

		for (k = 0; k < i; k++){
			d_x[i] -= d_m[i*n+k] * d_x[k];
		}
	}
}

	__global__ void LUPSolve2(double *d_m, int n, double *d_x) {
		int i, k;
		for (i = n - 1; i >= 0; i--) {
			for (k = i + 1; k < n; k++){
				d_x[i] -= d_m[i*n+k] * d_x[k];
			}

			d_x[i] /= d_m[i*n+i];
		}
	}

	/* INPUT: A - array of pointers to rows of a square matrix having dimension N
	*        Tol - small tolerance number to detect failure when the matrix is near degenerate
	* OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
	*        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
	*        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
	*        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
	*/
	__global__ void LUPDecompose(double *d_m, int n, double Tol, int *d_P) {
		int i, j, k, imax;
		double maxA, absA, temp;

		for (i = 0; i <= n; i++){
			d_P[i] = i; //Unit permutation matrix, P[N] initialized with N
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
				j = d_P[i];
				d_P[i] = d_P[imax];
				d_P[imax] = j;

				//pivoting rows of A
				for(j = 0; j < n; j++){
					temp = d_m[i*n+j];
					d_m[i*n + j] = d_m[imax*n+j];
					d_m[imax*n+j]=temp;
				}

				//counting pivots starting from n (for determinant)
				d_P[n]++;
			}

			for (j = i + 1; j < n; j++) {
				d_m[j*n+i] /= d_m[i*n+i];

				for (k = i + 1; k < n; k++){
					d_m[j*n+k] -= d_m[j*n+i] * d_m[i*n+k];
				}
			}
		}
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

		/*
		det = determinant(&m);

		if (det == 0 || (m.nrows != m.ncols)) {
		printf("ERROR: It is not possible to compute the inversion: determinant is equal to 0 or the matrix is not squared\n");
		fclose(mat);
		fclose(resultFile);
		free(m.mat);
		exit(1);
	}

	printf("Determinant of the matrix: %f \n", det);
	*/

	double *inverseM = (double*)malloc(m.ncols * m.nrows * sizeof(double));

	//cudaProfilerStart();

	/* Device variable, allocations, and transfers */
	double *d_m, *d_x, *d_b, *d_inverseM;
	int *d_P;

	cudaMalloc((void**)&d_m, m.nrows*m.ncols*sizeof(double));
	cudaMalloc((void**)&d_inverseM, m.nrows*m.ncols*sizeof(double));
	cudaMalloc((void**)&d_P, m.nrows**sizeof(int) + 1 * sizeof(int));
	cudaMalloc((void**)&d_x, m.nrows*sizeof(double));
	cudaMalloc((void**)&d_b, m.nrows*sizeof(double));

	cudaMemcpy(d_m, m.mat, m.nrows*m.ncols * sizeof(double), cudaMemcpyHostToDevice);

	int block_x = nelem / THREADS;
	if ((nelem) % THREADS != 0) {
		block_x++;
	}

	dim3 dimBlock(THREADS, 1, 1);
	dim3 dimGrid(block_x, 1, 1);

	dim3 dimBlock2(THREADS, THREADS);
	dim3 dimGrid2((m.ncols+dimBlock2.x-1)/dimBlock2.x, (m.nrows+dimBlock2.y-1)/dimBlock2.y);

	t = clock();

	LUPDecompose <<<(1,1,1), (1,1,1)>>>(d_m, m.nrows, 1E-3, d_P);
	LUPSolve <<<dimGrid, dimBlock>>>(d_m, d_P, d_b, m.nrows, d_x);
	LUPSolve2 <<<(1,1,1), (1,1,1)>>>(d_m, m.nrows, d_x);
	LUPInvert <<<dimGrid2, dimBlock2>>>(d_m, d_P, m.nrows, d_inverseM);

	cudaDeviceSynchronize();

	t = clock() - t;

	//printf("\nDeterminant: %lf\n", LUPDeterminant(m.mat, P, m.nrows));

	cudaMemcpy(inverseM, d_inverseM, m.nrows*m.ncols * sizeof(double), cudaMemcpyDeviceToHost);

	//cudaProfilerStop();

	resultFile = fopen("inverse.txt", "w");
	storeMatrix(inverseM, m.nrows, resultFile);

	printf("\nElapsed time: %lf seconds\n", ((double)t) / CLOCKS_PER_SEC);

	fclose(mat);
	fclose(resultFile);

	cudaFree(d_P);
	cudaFree(d_inverseM);
	cudaFree(d_x);
	cudaFree(d_b);
	cudaFree(d_m);

	free(inverseM);
	free(m.mat);

	return 0;
}
