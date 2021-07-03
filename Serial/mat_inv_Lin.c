/* 
* Same as in mat_inv.c but in this version we have tested the locality principle through the linearisation of the matrices.
* Unfortunately this code appeared to be slow compared with the other probably for the compilator optimization.
* Starting pseudocodes/algorithm: https://www-dimat.unipv.it/sangalli/numerical_methods_eng_sciences/1%20-%20Linear%20Systems%20(part%20I).pdf 
*/

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
void printMatrix(double *squared_matrix, int n, FILE* file);
double determinant(double *l, double *u, int n, int *perm);
void forwardSubstitution(double *l, double *p, double *y, int column, int n);
void backwardSubstitution(double *u, double *y, double *a_inv, int column, int n);
void pivoting(double *a, double *p, int n, int *perm);
void decomposition(double *l, double *u, int n);

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
void printMatrix(double *squared_matrix, int n, FILE* file) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fprintf(file, "%lf ", squared_matrix[i * n + j]);
		}
		fprintf(file, "\n");
	}
}

/* 
* Because LU decomposition is used, det M = det LU = det L * det U.
* L and U are triangular so the determinant is calculated as the product of the diagonal elements
*/
double determinant(double *l, double *u, int n, int *perm) {
	int i;
	double det = 1;

	for(i=0; i<n; i++){
		det *= l[i * n + i] * u[i * n + i];
	}
	
	return pow(-1,perm[0]) * det; //it is necessary to multiply the obtained the with (-1) raised to the number of permutation occurred in pivoting
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
void pivoting(double *a, double *p, int n, int *perm) { 
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
			perm[0]++;
    		}
	}
	free(temp);
}

/* Perf LU decomposition of matrix M to obtain matrices L (lower) and U (upper) used to resolve and equation system throud BW and FW to obtain the inverse */
void decomposition(double *l, double *u, int n) {
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
	if(argc != 2) { //Checking parameters: 1.mat_inv.exe 2.matrix.txt
		printf("Parameters error.\n");
		exit(1);
	}

	FILE *mat, *resultFile;
	clock_t t;
	struct matrix m;
	int i, perm=0;

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
	
	//printf("\nThe matrix you have inserted is %dx%d and has %d elements\nPlease wait until computation are done...\n\n", n,n,n*n);
	
	//Create pivoting and inverse matrices and matrices initialization
	double *a_inv = (double*)malloc(n * n * sizeof(double));
	double *p = (double*)malloc(n * n * sizeof(double));
	double *l = (double*)malloc(n * n * sizeof(double));
	double *a_p = (double*)malloc(n * n * sizeof(double));
	double *u = (double*)malloc(n * n * sizeof(double));
	double *y = (double*)malloc(n * sizeof(double));
    
	memset(a_inv, 0, n * n * sizeof(double));
	memset(p, 0, n * n * sizeof(double));
	memset(l, 0, n * n * sizeof(double));
	memcpy(a_p, m.mat, n * n * sizeof(double));
	
	for (i = 0; i < n; i++) {
        	p[i * n + i] = 1;
		l[i * n + i] = 1;
    	}
   
   	/* START */
	t = clock();
	
	pivoting(a_p, p, n, &perm);
	memcpy(u, a_p, n * n * sizeof(double));	//Fill u using a_p elements
	
    	decomposition(l, u, n);
	
	double det = determinant(l, u, n, &perm);
	printf("Determinant: %lf\n", det);
	
	if(det == 0.0){
		printf("ERROR: It is not possible to compute the inversion: determinant is equal to 0\n");
		fclose(mat);
		free(p);		
		free(l);
		free(u);
		free(a_p);
		free(y);
		free(a_inv);		
		free(m.mat);
		exit(1);
	}
	
	//Finding the inverse, result is stored into a_inv
	for (i = 0; i < n; i++) {
        	forwardSubstitution(l, p, y, i, n); 			// y is filled
        	backwardSubstitution(u, y, a_inv, i, n);		// a_inv is filled
    	}
	t = clock() - t;
	/* END */
	
	resultFile = fopen("inverse.txt", "w");
	printMatrix(a_inv, n, resultFile);

	printf("Elapsed time: %lf seconds\n", ((double)t) / CLOCKS_PER_SEC);

	fclose(mat);
	fclose(resultFile);
	free(p);		
	free(l);
	free(u);
	free(a_p);
	free(y);
	free(a_inv);		
	free(m.mat);

	return 0;
}
