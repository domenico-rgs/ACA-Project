#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
void matrixMul(struct matrix* m1, struct matrix* m2, struct matrix* m3) {
	int i, j, k;

	m3->nrows = m1->nrows;
	m3->ncols = m2->ncols;

	//dynamically allocates the m3-matrix to be used to save the result
	m3->mat = (double*)malloc(m3->nrows * m3->ncols * sizeof(double));
	memset(m3->mat, 0, m3->nrows * m3->ncols * sizeof(double)); //set the memory allocated to the m3 matrix to zero

	clock_t t = clock();
	//classical multiplication done with a summation
	for (i = 0; i < m1->nrows; i++) {
		for (j = 0; j < m2->ncols; j++) {
			for (k = 0; k < m1->ncols; k++) { 
				m3->mat[i * m3->ncols + j] += m1->mat[i*m1->ncols+k] * m2->mat[k*m2->ncols+j]; //"Linearisation" of matrices is used to achieve faster access due to the principle of access locality.
			}
		}
	}
	t = clock() - t;
	printf("time: %f seconds\n", ((double)t)/CLOCKS_PER_SEC);
}

int main(int argc, char* argv[]) {
	if(argc != 3){ //1- exe name, 2- mat1, 3- mat2
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

	//Multiplication is permitted if m1 is m x n and m2 is n x p.
	if(m1.ncols != m2.nrows){
		printf("It is not possible to do matrix multiplication. Check matrix number of rows and cols.\n");
		fclose(mat1);
		fclose(mat2);
		exit(1);
	}

	readMatrix(&m1, mat1);
	readMatrix(&m2, mat2);

	t = clock();
	matrixMul(&m1, &m2, &m3);
	t = clock() - t; //total time spent in matrixMul

	resultFile = fopen("result.txt", "w");
	printMatrix(&m3, resultFile);

	printf("Elapsed time: %f seconds\n", ((double)t)/CLOCKS_PER_SEC);

	fclose(mat1);
	fclose(mat2);
	fclose(resultFile);
	free(m1.mat);
	free(m2.mat);
	free(m3.mat);
	return 0;
}
