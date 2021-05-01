#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct matrix {
	int ncol;
	int nrow;
	int* mat;
};

void readMatrix(struct matrix* m, char* filename);
void printMatrice(struct matrix* m, char* filename);
void matrixMul(struct matrix* m1, struct matrix* m2, struct matrix* m3);

void readMatrix(struct matrix* m, char* filename) {
	FILE* file;
	int i, j;

	file = fopen(filename, "r");

	fscanf(file, "%d %d", &m->nrow, &m->ncol);

	m->mat = (int*)malloc(m->ncol * m->nrow * sizeof(int));

	for (i = 0; i < m->nrow; i++) {
		for (j = 0; j < m->ncol; j++) {
			fscanf(file, "%d", &m->mat[i * m->ncol + j]);
		}
	}
	fclose(file);
}

void printMatrix(struct matrix* m) {
	FILE* file;
	int i, j;

	file = fopen("result.txt", "w");

	for (i = 0; i < m->nrow; i++) {
		for (j = 0; j < m->ncol; j++) {
			fprintf(file, "%d ", m->mat[i * m->ncol + j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void matrixMul(struct matrix* m1, struct matrix* m2, struct matrix* m3) {
	int i, j, k;

	m3->nrow = m1->nrow;
	m3->ncol = m2->ncol;

	m3->mat = (int*)malloc(m3->nrow * m3->ncol * sizeof(int));

	memset(m3->mat, 0, m3->nrow * m3->ncol * sizeof(int));

	for (i = 0; i < m1->nrow; i++) {
		for (j = 0; j < m2->ncol; j++) {
			for (k = 0; k < m1->ncol; k++) {
				m3->mat[i * m3->ncol + j] += m1->mat[i*m1->ncol+k] * m2->mat[k*m2->ncol+j];
			}
		}
	}
}

int main(int argc, char* argv[]) {
  clock_t t;
	struct matrix m1, m2, m3;

	readMatrix(&m1, argv[1]);
	readMatrix(&m2, argv[2]);

  t = clock();
	matrixMul(&m1, &m2, &m3);
  t = clock() - t;
	printMatrix(&m3);

  printf("Elapsed time: %f seconds", ((double)t)/CLOCKS_PER_SEC);

	free(m1.mat);
	free(m2.mat);
	free(m3.mat);
	return 0;
}
