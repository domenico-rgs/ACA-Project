#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cuda_runtime.h"
//#include "cuda_profiler_api.h"

#define THREADS 32 //In each block THREADS*THREADS threads

struct matrix {
	int ncols;
	int nrows;
	double* mat;
};

void readMatrix(struct matrix* m, FILE* file);
void printMatrix(struct matrix* m, FILE* file);
__global__ void matrixMul(double *d_m1, double *d_m2, double *d_m3, int row1, int row2, int col1, int col2);

/*
 * The file which contains a matrix has in its first row the dimensions
 * then using fscanf each element of the matrix is stored on the memory allocated dynamically
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

/* The opposite operation of readMatrix. Stores a matrix into a file, element by element */
void printMatrix(struct matrix* m, FILE* file) {
	int i, j;

	for (i = 0; i < m->nrows; i++) {
		for (j = 0; j < m->ncols; j++) {
			fprintf(file, "%lf ", m->mat[i * m->ncols + j]);
		}
		fprintf(file, "\n");
	}
}

/*
 * Performs the multiplication operation between the matrices m1 and m2.
 * The result will be stored in the matrix m3.
 * The algorithm is practically the one that can be found here: https://en.wikipedia.org/wiki/Matrix_multiplication#Definition
 * The grid is sized in a 2D way as each blocks whose dimension is THREADSxTHREADS, this helps in associating to a single thread and element of the result matrix to which calcutate a summation exploiting the numerous core of a GPU.
*/
__global__ void matrixMul(double *d_m1, double *d_m2, double *d_m3, int row1, int row2, int col1, int col2){
	//coordinates of each thread inside the grind and so inside a block
	int i = blockIdx.y*blockDim.y+threadIdx.y;
	int j = blockIdx.x*blockDim.x+threadIdx.x;

	double sum = 0;
	int k;

	//the two previous for cycle are substituted by the matrix of threads
	if ((i < row1) && (j < col2)){
		for(k = 0; k<col1; k++){
			sum += d_m1[i*col1+k]*d_m2[k*col2+j];
		}
		d_m3[i*col2+j]=sum;
  }
}

int main(int argc, char* argv[]) {
	if(argc != 3){ //1- exe name, 2- mat1.txt, 3- mat2.txt
		printf("Parameter error.\n");
		exit(1);
	}

	FILE *mat1, *mat2, *resultFile;
	clock_t t;
	struct matrix m1, m2, m3;

	mat1 = fopen(argv[1], "r");
	mat2 = fopen(argv[2], "r");
	fscanf(mat1, "%d %d", &m1.nrows, &m1.ncols);
	fscanf(mat2, "%d %d", &m2.nrows, &m2.ncols);

	/* Multiplication is permitted if m1 is m x n and m2 is n x p, m1 must have the same number of column of the rows of m2 matrix */
	if(m1.ncols != m2.nrows){
		printf("It is not possible to do matrix multiplication. Check matrices number of rows and cols.\n");
		fclose(mat1);
		fclose(mat2);
		exit(1);
	}

	readMatrix(&m1, mat1);
	readMatrix(&m2, mat2);

	//M3 initilization
	m3.nrows=m1.nrows;
	m3.ncols=m2.ncols;
	m3.mat = (double*)malloc(m3.ncols * m3.nrows * sizeof(double));

	//cudaProfilerStart();

	/* Device variable, allocations, and transfers */
	double *d_m1, *d_m2, *d_m3;
	cudaMalloc((void**)&d_m1, m1.nrows*m1.ncols*sizeof(double));
	cudaMalloc((void**)&d_m2, m2.nrows*m2.ncols*sizeof(double));
	cudaMalloc((void**)&d_m3, m3.nrows*m3.ncols*sizeof(double));

	cudaMemcpy(d_m1, m1.mat, m1.nrows*m1.ncols * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m2, m2.mat, m2.nrows*m2.ncols * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemset(d_m3, 0, m3.nrows*m3.ncols*sizeof(double));

	dim3 dimBlock(THREADS, THREADS); //each block is THREADxTHREAD
	dim3 dimGrid((m2.ncols+dimBlock.x-1)/dimBlock.x, (m1.nrows+dimBlock.y-1)/dimBlock.y); //the grid is allocated so the necessary amount of blocks are initialized

	t = clock();
	matrixMul <<<dimGrid, dimBlock>>>(d_m1, d_m2, d_m3, m1.nrows, m2.nrows, m1.ncols, m2.ncols);
	cudaDeviceSynchronize();
	t = clock() - t; //total time spent in matrixMul (wall clock time)

	cudaMemcpy(m3.mat, d_m3, m3.nrows*m3.ncols * sizeof(double), cudaMemcpyDeviceToHost);

	//cudaProfilerStop();

	resultFile = fopen("result.txt", "w");
	printMatrix(&m3, resultFile);

	printf("Elapsed time: %.5lf seconds\n", ((double)t)/CLOCKS_PER_SEC);

	fclose(mat1);
	fclose(mat2);
	fclose(resultFile);
	cudaFree(d_m1);
	cudaFree(d_m2);
	cudaFree(d_m3);
	free(m1.mat);
	free(m2.mat);
	free(m3.mat);

	return 0;
}
