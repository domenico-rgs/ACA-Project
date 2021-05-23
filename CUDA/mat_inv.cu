#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cuda_runtime.h"
//#include "cuda_profiler_api.h"

#define THREADS 32

struct matrix {
	int ncols;
	int nrows;
	double* mat;
};

void readMatrix(struct matrix* m, FILE* file);
void printMatrix(double *squared_matrix, int n, FILE* file);
__global__ void determinant(double *l, double *u, int n, double *det);
void forwardSubstitution(double *l, double *p, double *y, int column, int n);
void backwardSubstitution(double *u, double *y, double *a_inv, int column, int n);
void pivoting(double *a, double *p, int n);
void lu(double *l, double *u, int n);
__device__ double atomicMul(double* address, double val);

/* Reads a matrix from a file and stores it into the appropriate structure. */
void readMatrix(struct matrix* m, FILE* file) {
	int i, j;

	m->mat = (double*)malloc(m->ncols * m->nrows * sizeof(double));

	for (i = 0; i < m->nrows; i++) {
		for (j = 0; j < m->ncols; j++) {
			fscanf(file, "%lf", &m->mat[i * m->ncols + j]);
		}
	}
}

/* Stores a matrix into the file passed as argument */
void printMatrix(double *squared_matrix, int n, FILE* file) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fprintf(file, "%lf ", squared_matrix[i * n + j]);
		}
		fprintf(file, "\n");
	}
}

#if __CUDA_ARCH__ < 600
__device__ double atomicMul(double* address, double val){
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val *
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

/* Becaute LU decomposition is used  det M = det LU = det L * det U, L and U are triangular
   so the determinant is calculated as the product of the diagonal elements
 */
__global__ void determinant(double *l, double *u, int n, double *d_det) {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	double det = 1;

	if(i==1){
		d_det[0]=1;
	}
	__syncthreads();

	if(i<n){
			det *= l[i * n + i] * u[i * n + i];
	}

	atomicMul(d_det, det);
}

/* Since L is a lower triangular matrix forward substitution is used to perform the calculus of Lx=y */
void forwardSubstitution(double *l, double *p, double *y, int column, int n) {
	int i, j;
	double sum = 0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            sum = sum + l[i * n + j] * y[j];
        }
        y[i] = (p[i * n + column] - sum) / l[i * n + i];
        sum = 0;
    }
}

/* Since U is an upper triangular matrix backward substitution is used to perform the calculus of Ux=y */
void backwardSubstitution(double *u, double *y, double *a_inv, int column, int n) {
    int i, j;
	double sum;

    a_inv[(n-1)*n+column] = y[n-1] / u[(n-1) * n + (n-1)];
    for (i = n - 2; i >= 0; i--) {
		sum = y[i];
        for (j = n - 1; j > i; j--) {
           sum = sum - u[i * n + j] * a_inv[j*n+column];
        }
        a_inv[i*n+column] = sum / u[i * n + i];
        sum = 0;
    }
}

/* Even if det(M)!=0, pivoting is performed to be sure that L and U are correctly upper and lower triangular matrix */
void pivoting(double *a, double *p, int n) {
    int j, k;
	int isMaximum = 0;
    double *temp = (double*)malloc(n * sizeof(double));

    // k is column and j is row
	for (k = 0; k < n-1; k++) {
    	int imax = k;
        for (j = k; j < n; j++) {
			if (a[j * n + k] > a[imax * n + k]) {  // finding the maximum index
				imax = j;
                isMaximum = 1;
            }
        }
        if (isMaximum == 1) {
        	// swapping a[k] and a[imax]
			memcpy(temp, &a[k*n], n * sizeof(double));
			memcpy(&a[k*n], &a[imax*n], n * sizeof(double));
			memcpy(&a[imax*n], temp, n * sizeof(double));

			// swapping p[k] and p[imax]
			memcpy(temp, &p[k*n], n * sizeof(double));
			memcpy(&p[k*n], &p[imax*n], n * sizeof(double));
			memcpy(&p[imax*n], temp, n * sizeof(double));

        	isMaximum = 0;
    	}
	}
	free(temp);
}

/* Perf LU decomposition of matrix M*/
void lu(double *l, double *u, int n) {
    int i, j, k;

	for (k = 0; k < n; k++) {
        for (i = k + 1; i < n; i++) {
            l[i * n + k] = u[i * n + k] / u[k * n + k];
            for (j = k; j < n; j++) {
                u[i * n + j] = u[i * n + j] - l[i * n + k] * u[k * n + j];
            }
        }
    }
}

int main(int argc, char* argv[]) {
	if(argc != 2) { //Checking parameters: 1.mat_inv.exe 2.matrix
		printf("Parameters error.\n");
		exit(1);
	}

	printf("This program compute the inverse of a squared matrix using only one thread\nPlease wait until computation are done...\n");

	FILE *mat, *resultFile;
	clock_t t;
	struct matrix m;
	int i;

	mat = fopen(argv[1], "r");
	fscanf(mat, "%d %d", &m.nrows, &m.ncols);
	readMatrix(&m, mat);

	if (m.nrows != m.ncols) {
		printf("ERROR: It is not possible to compute the inversion: the matrix is not squared\n");
		fclose(mat);
		free(m.mat);
		exit(1);
	}

	int n = m.nrows; //matrix order (m is squared)

	//Create pivoting and inverse matrices
	double *a_inv = (double*)malloc(n * n * sizeof(double));
	double *p = (double*)malloc(n * n * sizeof(double));
	double *l = (double*)malloc(n * n * sizeof(double));
	double *a_p = (double*)malloc(n * n * sizeof(double));
	double *u = (double*)malloc(n * n * sizeof(double));
	double *y = (double*)malloc(n * sizeof(double));

	//Matrices initialization
	memset(a_inv, 0, n * n * sizeof(double));
	memset(p, 0, n * n * sizeof(double));
	memset(l, 0, n * n * sizeof(double));
	memcpy(a_p, m.mat, n * n * sizeof(double));

	for (i = 0; i < n; i++) {
        p[i * n + i] = 1;
		l[i * n + i] = 1;
    }

	/* Device variable, allocations, and transfers */
	double *d_l, *d_u, *d_det;
	double det;

	cudaMalloc((void**)&d_l, n*n*sizeof(double));
	cudaMalloc((void**)&d_u, n*n*sizeof(double));
   	cudaMalloc((void**)&d_det, sizeof(double));

	int block_x = n / THREADS;
	if ((n) % THREADS != 0) {
		block_x++;
	}

	dim3 dimBlock(THREADS, 1, 1);
	dim3 dimGrid(block_x, 1, 1);

	t = clock();
	pivoting(a_p, p, n);

	memcpy(u, a_p, n * n * sizeof(double));	//Fill u using a_p elements

    lu(l, u, n);

	cudaMemcpy(d_l, l, n * n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_u, u, n * n * sizeof(double), cudaMemcpyHostToDevice);

	determinant <<<dimGrid, dimBlock>>>(d_l, d_u, n, d_det);

	cudaMemcpy(&det, d_det, sizeof(double), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	printf("Determinant: %lf\n", det);
	if(det == 0.0){
		printf("ERROR: It is not possible to compute the inversion: the matrix is not squared\n");
		fclose(mat);
		cudaFree(d_l);
		cudaFree(d_u);
		cudaFree(d_det);
		free(p);
		free(l);
		free(u);
		free(a_p);
		free(y);
		free(a_inv);
		free(m.mat);
		exit(1);
	}

	/* Finding the inverse, result is stored into a_inv */
	for (i = 0; i < n; i++) {
        forwardSubstitution(l, p, y, i, n); 			// y is filled
        backwardSubstitution(u, y, a_inv, i, n);		// a_inv is filled
    }
	t = clock() - t;

	//cudaProfilerStop();

	resultFile = fopen("inverse.txt", "w");
	printMatrix(a_inv, n, resultFile);

	printf("\nElapsed time: %lf seconds\n", ((double)t) / CLOCKS_PER_SEC);

	fclose(mat);
	fclose(resultFile);

	cudaFree(d_l);
	cudaFree(d_u);
	cudaFree(d_det);


	free(p);
	free(l);
	free(u);
	free(a_p);
	free(y);
	free(a_inv);
	free(m.mat);

	return 0;
}
