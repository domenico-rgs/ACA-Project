#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cuda_runtime.h"
//#include "cuda_profiler_api.h"

#define TILE_WIDTH 2

struct matrix {
	int ncols;
	int nrows;
	double* mat;
};

void readMatrix(struct matrix* m, FILE* file);
void printMatrice(struct matrix* m, FILE* file);
void matrixMul(struct matrix* m1, struct matrix* m2, struct matrix* m3);

/*
Knowing the number of rows and columns,
it reads a matrix from a file and stores it in the appropriate structure.
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
The opposite operation to readMatrix. Saves a matrix in the file given as argument
*/
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
Performs the multiplication operation between the matrices m1 and m2.
The result will be stored in the matrix m3.
*/
__global__ void matrixMul(double *d_m1, double *d_m2, double *d_m3, int width){
  int row = blockIdx.y*TILE_WIDTH + threadIdx.y;
  int col = blockIdx.x*TILE_WIDTH + threadIdx.x;

  double Pvalue=0;
  for(int k=0; k<width; k++){
    Pvalue += d_m1[row*width+k]+d_m2[k*width+col];
  }
  d_m3[row*width+col]=Pvalue;
}

int main(int argc, char* argv[]) {
	if(argc != 3){ //1- exe name, 2- mat1, 3- mat2
		printf("Parameter error.");
		exit(1);
	}

	FILE *mat1, *mat2, *resultFile;
	clock_t t;
	struct matrix m1, m2, m3;

	mat1 = fopen(argv[1], "r");
	mat2 = fopen(argv[2], "r");
	fscanf(mat1, "%d %d", &m1.nrows, &m1.ncols);
	fscanf(mat2, "%d %d", &m2.nrows, &m2.ncols);

	//Multiplication is permitted if m1 is m x n and m2 is n x p.
	if(m1.ncols != m2.nrows){
		printf("It is not possible to do matrix multiplication. Check matrix number of rows and cols.");
		fclose(mat1);
		fclose(mat2);
		exit(1);
	}

	readMatrix(&m1, mat1);
	readMatrix(&m2, mat2);

  //cudaProfilerStart();

  /* Device variable, allocations, and transfers */
  double *d_m1, *d_m2, *d_m3;
  cudaMalloc((void**)&d_m1, sizeof(double) * m1.nrows*m1.ncols);
  cudaMalloc((void**)&d_m2, sizeof(double) * m2.nrows*m2.ncols);
  cudaMalloc((void**)&d_m3, sizeof(double) * m3.nrows*m3.ncols);

  cudaMemcpy(d_m1, m1.mat, m1.nrows*m1.ncols * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_m2, m2.mat, m2.nrows*m2.ncols * sizeof(double), cudaMemcpyHostToDevice);

  dim3 dimGrid(m1.nrows/TILE_WIDTH, m1.nrows/TILE_WIDTH);
  dim3 dimBlock(TILE_WIDTH, TILE_WIDTH);

	t = clock();
  matrixMul <<<dimGrid, dimBlock>>>(d_m1, d_m2, d_m3, m1.nrows);
	t = clock() - t; //total time spent in matrixMul

	resultFile = fopen("result.txt", "w");
  cudaMemcpy(d_m3, m3.mat, m3.nrows*m3.ncols * sizeof(double), cudaMemcpyDeviceToHost);
	printMatrix(&m3, resultFile);

	printf("Elapsed time: %f seconds", ((double)t)/CLOCKS_PER_SEC);

	fclose(mat1);
	fclose(mat2);
	fclose(resultFile);

  cudaFree(d_m1);
  cudaFree(d_m2);
  cudaFree(d_m3);

	free(m1.mat);
	free(m2.mat);
	free(m3.mat);

  //cudaProfilerStop();
	return 0;
}
