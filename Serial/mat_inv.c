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
void storeMatrix(double *squared_matrix, int n, FILE* file);
double findDeterminant(double *a, int *p, int n);
void forwardSubstitution(double *l, double *p, double *y, int column, int n);
void backwardSubstitution(double *u, double *y, double *a_inv, int column, int n);
void pivoting(double *a, double *p, int n);
void lu(double *l, double *u, int n);
void computeInverse(double *a, double *a_inv, int n);

/*
 *	Reads a matrix from a file and stores it into the appropriate structure.
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

/* TODO */
double findDeterminant(double *a, int *p, int n) {
	int i;
	double det = a[0];

	for (i = 1; i < n; i++) {
		det *= a[i*n + i];
	}

	return (p[n] - n) % 2 == 0 ? det : -det;
}

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

void pivoting(double *a, double *p, int n) { 
    int j, k;
	int isMaximum = 0; 
    double *temp = (double*)malloc(n * sizeof(double));	// to perform memory swappings
    
    // k is column and j is row
	for (k = 0; k < n-1; k++) {   
    	int imax = k;
        for (j = k; j < n; j++) { 	
            // finding the maximum index
			if (a[j * n + k] > a[imax * n + k]) { 
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

void computeInverse(double *a, double *a_inv, int n) {
    int i, j;
    
    /* Memory allocation of temporary matrix */
    double *p = (double*)malloc(n * n * sizeof(double));
	double *l = (double*)malloc(n * n * sizeof(double));
	double *a_p = (double*)malloc(n * n * sizeof(double));
	double *u = (double*)malloc(n * n * sizeof(double));
	double *y = (double*)malloc(n * sizeof(double));
    
	/* Creating pivoting matrix and setting inverse matrix to 0 */
	for (i = 0; i < n; i++) {
        p[i * n + i] = 1;
        for (j = 0; j < n; j++) {
            a_inv[i*n+j] = 0;
            if (i != j) {
                p[i * n + j] = 0;
            }
        }
    }
    
    /* Creating matrix l and copying a into a_p */
    for (i = 0; i < n; i++) {
        l[i * n + i] = 1;
        for (j = 0; j < n; j++) {
            a_p[i * n + j] = a[i*n + j];
            if (i != j)
                l[i * n + j] = 0;
        }
    }
    
    pivoting(a_p, p, n);

	/* Creating matrix u using a_p */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            u[i * n + j] = a_p[i * n + j];
        }
    }

	/* Performing LU decomposition */
    lu(l, u, n);
    
    /* Finding the inverse, result is stored into a_inv */
	for (i = 0; i < n; i++) {
        forwardSubstitution(l, p, y, i, n); 			// y is filled
        backwardSubstitution(u, y, a_inv, i, n);		// a_inv is filled
    }
    
    free(p);		
	free(l);
	free(u);
	free(a_p);
	free(y);
}

int main(int argc, char* argv[]) {
	/* Checking parameters: 1.mat_inv.exe 2.matrix */
	if(argc != 2) {
		printf("Parameters error.\n");
		exit(1);
	}

	FILE *mat, *resultFile;
	clock_t t;
	struct matrix m;
	//int i, j;

	mat = fopen(argv[1], "r");
	fscanf(mat, "%d %d", &m.nrows, &m.ncols);
	readMatrix(&m, mat);

	if (m.nrows != m.ncols) {
		printf("ERROR: It is not possible to compute the inversion: the matrix is not squared\n");
		fclose(mat);
		free(m.mat);
		exit(1);
	}
	
	/*	
	double det = findDeterminant(m.mat, p, m.nrows);
	printf("\nDeterminant: %lf\n", det);
	if (det == 0.0) {
		printf("ERROR: It is not possible to compute the inversion\n");
		fclose(mat);
		free(*a1);		
		free(m.mat);
		exit(1);
	}
	*/
	
	double *a_inv = (double*)malloc(m.nrows * m.ncols * sizeof(double));

    
	t = clock();
	computeInverse(m.mat, a_inv, m.nrows);
	t = clock() - t;
	
	resultFile = fopen("inverse.txt", "w");
	storeMatrix(a_inv, m.nrows, resultFile);

	printf("\nElapsed time: %lf seconds\n", ((double)t) / CLOCKS_PER_SEC);

	fclose(mat);
	fclose(resultFile);
	free(a_inv);		
	free(m.mat);

	return 0;
}