#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct matrix {
	int ncols;
	int nrows;
	int* mat;
};

void readMatrix(struct matrix* m, FILE* file);
void printMatrice(struct matrix* m, FILE* file);
void matrixMul(struct matrix* m1, struct matrix* m2, struct matrix* m3);

void readMatrix(struct matrix* m, FILE* file) {
	int i, j;

	m->mat = (int*)malloc(m->ncols * m->nrows * sizeof(int));

	for (i = 0; i < m->nrows; i++) {
		for (j = 0; j < m->ncols; j++) {
			fscanf(file, "%d", &m->mat[i * m->ncols + j]);
		}
	}
}

void printMatrix(struct matrix* m, FILE* file) {
	int i, j;

	for (i = 0; i < m->nrows; i++) {
		for (j = 0; j < m->ncols; j++) {
			fprintf(file, "%d ", m->mat[i * m->ncols + j]);
		}
		fprintf(file, "\n");
	}
}

void matrixMul(struct matrix* m1, struct matrix* m2, struct matrix* m3) {
	int i, j, k;

	m3->nrows = m1->nrows;
	m3->ncols = m2->ncols;

	m3->mat = (int*)malloc(m3->nrows * m3->ncols * sizeof(int));

	memset(m3->mat, 0, m3->nrows * m3->ncols * sizeof(int));

	for (i = 0; i < m1->nrows; i++) {
		for (j = 0; j < m2->ncols; j++) {
			for (k = 0; k < m1->ncols; k++) {
				m3->mat[i * m3->ncols + j] += m1->mat[i*m1->ncols+k] * m2->mat[k*m2->ncols+j];
			}
		}
	}
}

int main(int argc, char* argv[]) {
	if(argc != 3){
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

	if(m1.ncols != m2.nrows){
		printf("It is not possible to do matrix multiplication. Check matrix number of rows and cols.");
		fclose(mat1);
		fclose(mat2);
		exit(1);
	}
	readMatrix(&m1, mat1);
	readMatrix(&m2, mat2);

  t = clock();
	matrixMul(&m1, &m2, &m3);
  t = clock() - t;

	resultFile = fopen("result.txt", "w");
	printMatrix(&m3, resultFile);

  printf("Elapsed time: %f seconds", ((double)t)/CLOCKS_PER_SEC);

	fclose(mat1);
	fclose(mat2);
	fclose(resultFile);
	free(m1.mat);
	free(m2.mat);
	free(m3.mat);
	return 0;
}
