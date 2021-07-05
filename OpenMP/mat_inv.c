#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

void readMatrix(double** matrix, FILE* file, int n);
void printMatrix(double **matrix, int n, FILE* file);
double determinant(double **l, double **u, int n, int *perm);
void forwardSubstitution(double **l, double **p, double *y, int column, int n);
void backwardSubstitution(double **u, double *y, double **a_inv, int column, int n);
void pivoting(double **a, double **p, int n, int *perm);
void decomposition(double **l, double **u, int n);

/*
 * The file which contains a matrix has in its first row the dimensions
 * then using fscanf each element of the matrix is stored on the memory allocated dynamically
*/
void readMatrix(double** matrix, FILE* file, int n) {
	int i, j;

	for (i = 0; i < n; i++){
		matrix[i] = (double*)malloc(n * sizeof(double));
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fscanf(file, "%lf", &matrix[i][j]);
		}
	}
}

/* The opposite operation of readMatrix. Stores a matrix into a file, element by element */
void printMatrix(double **matrix, int n, FILE* file) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fprintf(file, "%lf ", matrix[i][j]);
		}
		fprintf(file, "\n");
	}
}

/* 
* Because LU decomposition is used, det M = det LU = det L * det U.
* L and U are triangular so the determinant is calculated as the product of the diagonal elements
*/
double determinant(double **l, double **u, int n, int *perm) {
	int i;
	double det = 1;
	
	#pragma omp parallel
	{
		#pragma omp for reduction(*: det) //to speedup computation we apply the reduction to det. In this way the iteration are distributed among threads and at the end the final det value is calculated with a parallel reduction
		for(i = 0; i < n; i++) 
			det *= l[i][i] * u[i][i];
	}
	
	return pow(-1, perm[0]) * det; //it is necessary to multiply the obtained the with (-1) raised to the number of permutation occurred in pivoting
}

/* Since L is a lower triangular matrix forward substitution is used to perform the calculus of Lx=y */
void forwardSubstitution(double **l, double **p, double *y, int column, int n) {
	int i, j;
	double sum = 0;
	
	for (i = 0; i < n; i++) {
		//#pragma omp parallel for reduction(+:sum) //another way of parallelizing spreading the work for each linear equation among threads but it is better to do a parallelizazione as done in the for loop in main
		for (j = 0; j < i; j++) {
			sum = sum + l[i][j] * y[j];
        	}
        	y[i] = (p[i][column] - sum) / l[i][i];
        	sum = 0;
	}
}

/* Since U is an upper triangular matrix backward substitution is used to perform the calculus of Ux=y */
void backwardSubstitution(double **u, double *y, double **a_inv, int column, int n) {
	int i, j;
	double sum;
    
	a_inv[n-1][column] = y[n-1] / u[n-1][n-1];
	for (i = n - 2; i >= 0; i--) {
		sum = y[i];
		//#pragma omp parallel for reduction(+:sum) //as above
		for (j = n - 1; j > i; j--) {
			sum = sum - u[i][j] * a_inv[j][column];
       		 }
        	a_inv[i][column] = sum / u[i][i];
       		sum = 0;
  	 }
}

/* Even if det(M)!=0, pivoting is performed to be sure that L and U are correctly upper and lower triangular matrix */
void pivoting(double **a, double **p, int n, int *perm) {
	int j, k;
	int isMaximum = 0; 
	double *temp = (double*)malloc(n * sizeof(double));
    	
	// k is column and j is row
	for (k = 0; k < n-1; k++) {   
		int imax = k;
    		for (j = k; j < n; j++) { 	
			if (a[j][k] > a[imax][k]) {  // finding the maximum index
				imax = j;
	            		isMaximum = 1;
	       		 }
    		}
    		if (isMaximum == 1) {
    			// swapping a[k] and a[imax]
			memcpy(temp, a[k], n * sizeof(double));
			memcpy(a[k], a[imax], n * sizeof(double));
			memcpy(a[imax], temp, n * sizeof(double));
			
			// swapping p[k] and p[imax]
			memcpy(temp, p[k], n * sizeof(double));
			memcpy(p[k], p[imax], n * sizeof(double));
			memcpy(p[imax], temp, n * sizeof(double));
			
	    		isMaximum = 0;
	    		
			perm[0]++;
		}
	}
	free(temp);
}

/* Perf LU decomposition of matrix M to obtain matrices L (lower) and U (upper) used to resolve and equation system throud BW and FW to obtain the inverse */
void decomposition(double **l, double **u, int n) {
    int i, j, k;
    
    	#pragma omp parallel private (i,j,k)
    	{
		for (k = 0; k < n; k++) { //u is shared, parallelizing this for loop is not possibile because otherwise there will be data races with u
			#pragma omp for schedule(static) //computation of rows of l and u are done in parallel with different threads
			for (i = k + 1; i < n; i++) {
            			l[i][k] = u[i][k] / u[k][k];
            			for (j = k; j < n; j++) {
            				u[i][j] = u[i][j] - l[i][k] * u[k][j];
           			 }	
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
	double t;
	int i, rows, cols, perm = 0;

	mat = fopen(argv[1], "r");
	fscanf(mat, "%d %d", &rows, &cols);
	
	if (rows != cols) {
		printf("ERROR: It is not possible to compute the inversion: the matrix is not squared\n");
		fclose(mat);
		exit(1);
	}
	
	int n = rows; //matrix order (m is squared)
	double **m = (double **)malloc(n * sizeof(double*));
	readMatrix(m, mat, n);

	//printf("\nThe matrix you have inserted is %dx%d and has %d elements\nPlease wait until computation are done...\n\n", n, n, n * n);
	
	/* Create pivoting and inverse matrices and matrices initialization */
	double **a_inv = (double **)malloc(n * sizeof(double*));
	double **p = (double **)malloc(n * sizeof(double*));
	double **l = (double **)malloc(n * sizeof(double*));
	double **a_p = (double **)malloc(n * sizeof(double*));
	double **u = (double **)malloc(n * sizeof(double*));

	#pragma omp parallel for //spread initialization work but it not speedup a lot the code so it was ignored in wall clock time measures
	for(i = 0; i < n; i++) { 
		a_inv[i] = (double *)malloc(n * sizeof(double));
		p[i] = (double *)malloc(n * sizeof(double));
		l[i] = (double *)malloc(n * sizeof(double));
		a_p[i] = (double *)malloc(n * sizeof(double));
		u[i] = (double *)malloc(n * sizeof(double));
		
		memset(a_inv[i], 0, n * sizeof(double));
		memset(p[i], 0, n * sizeof(double));
		memset(l[i], 0, n * sizeof(double));
		memset(u[i], 0, n * sizeof(double));
		memcpy(a_p[i], m[i], n * sizeof(double));
		
		p[i][i] = 1;
		l[i][i] = 1;
	}
    
   	/* START */
	t = omp_get_wtime();	
	pivoting(a_p, p, n, &perm);

	#pragma omp parallel for
	for (i = 0; i < n; i++){
		memcpy(u[i], a_p[i], n * sizeof(double));	// Fill u using a_p elements
	}

   	 decomposition(l, u, n);

	double det = determinant(l, u, n, &perm);
	printf("Determinant: %lf\n", det);
	
	if(det == 0.0) {
		printf("ERROR: It is not possible to compute the inversion: determinant is equal to 0\n");
		fclose(mat);
		for (i = 0; i < n; i++) {
			free(p[i]);
			free(a_p[i]);
			free(u[i]);
			free(l[i]);
			free(a_inv[i]);
			free(m[i]);
		}
		free(p);		
		free(l);
		free(u);
		free(a_p);
		free(a_inv);		
		free(m);
		exit(1);
	}

	/* Finding the inverse, result is stored into a_inv */
	#pragma omp parallel shared(a_inv) private(i)
	{
		#pragma omp for schedule(dynamic)
		for (i = 0; i < n; i++) {
			double *y = (double*)malloc(n * sizeof(double));
	        	forwardSubstitution(l, p, y, i, n); 			// y is filled
	        	backwardSubstitution(u, y, a_inv, i, n);		// a_inv is filled
	        	free(y);
	    	}
	}
	t = omp_get_wtime() - t;
	/* STOP */
	resultFile = fopen("inverse.txt", "w");
	printMatrix(a_inv, n, resultFile);

	printf("Elapsed time: %lf seconds\n", t);

	fclose(mat);
	fclose(resultFile);
	
	for (i = 0; i < n; i++) {
		free(p[i]);
		free(a_p[i]);
		free(u[i]);
		free(l[i]);
		free(a_inv[i]);
		free(m[i]);
	}
	free(p);		
	free(l);
	free(u);
	free(a_p);
	free(a_inv);		
	free(m);

	return 0;
}
