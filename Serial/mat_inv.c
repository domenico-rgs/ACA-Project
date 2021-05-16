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
void storeMatrix(struct matrix* m, FILE* file);
void matrixInversion(struct matrix* m1, struct matrix* m2);
double computeDeterminant(struct matrix* m);

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
void storeMatrix(struct matrix* m, FILE* file) {
	int i, j;

	for (i = 0; i < m->nrows; i++) {
		for (j = 0; j < m->ncols; j++) {
			fprintf(file, "%lf ", m->mat[i * m->ncols + j]);
		}
		fprintf(file, "\n");
	}
}

/* https://www.codesansar.com/c-programming-examples/matrix-determinant.htm */
double computeDeterminant(struct matrix* m) {
	int i, j, k, ratio;
	double det = 1;
	double* mat;

	mat = (double*)malloc(m->ncols * m->nrows * sizeof(double));
	memcpy(mat, m->mat, m->nrows * m->ncols * sizeof(double)); //necessary in order not to change the read matrix into an upper triangular one

	/* Using Gauss Elimination technique for transforming matrix to upper triangular matrix */
	for (i = 0; i < m->nrows; i++) {
	 	if (mat[i * m->ncols + i] == 0) {
			free(mat);
			return HUGE_VAL;
		}
		
		for (j = i+1 ; j < m->nrows; j++) {
			ratio = mat[j * m->ncols + i] / mat[i * m->ncols + i];

			for (k = 0; k < m->nrows; k++) {
				mat[j * m->ncols + k] = mat[j * m->ncols + k] - ratio * mat[i * m->ncols + k];
			}
		}
	}

	printf("\nUpper Triangular Matrix: \n");
	 for(i=0;i< m->nrows;i++)
	 {
		  for(j=0;j< m->nrows;j++)
		  {
			   printf("%0.2f\t",mat[i*m->ncols+j]);
		  }
		  printf("\n");
	 }

	/* Finding determinant by multiplying elements in principal diagonal elements */
	for (i = 0; i < m->nrows; i++) {
        det = det * mat[i * m->ncols + i];
    }

    free(mat);
    return det;
}

/*
	Performs the inversion of the matrix m by using Jacobi algorithm
*/
void matrixInversionJacobi(struct matrix* m) {
	int i, j, k;

	/* Jacobi Algorithm */

}

/*
	Performs the inversion of the matrix m by using Gauss-Seidel algorithm
*/
void matrixInversionGauss(struct matrix* m) {
	int i, j, k;

	/* Gauss-Seidel Algorithm */
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

	det = computeDeterminant(&m);

	/* Checking if it is possible to perform the matrix inversion */
	if(det == HUGE_VAL){
		printf("ERROR: Mathematical error in determinant calculation\n");
		fclose(mat);
		fclose(resultFile);
		free(m.mat);
		exit(1);
	}else if (det == 0 || (m.nrows != m.ncols)) {
		printf("ERROR: It is not possible to compute the inversion: determinant is equal to 0 or the matrix is not squared\n");
		fclose(mat);
		fclose(resultFile);
		free(m.mat);
		exit(1);
	}

	printf("Determinant of the matrix: %f \n", det);

	t = clock();
	matrixInversionJacobi(&m);
	t = clock() - t;

	resultFile = fopen("inverse.txt", "w");
	storeMatrix(&m, resultFile);

	printf("\nElapsed time: %lf seconds\n", ((double)t) / CLOCKS_PER_SEC);

	fclose(mat);
	fclose(resultFile);
	free(m.mat);

	return 0;
}
