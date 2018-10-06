#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl_lapacke.h"

double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s\n", matrix);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if (bref[i]!=b[i]) return 0;
    }
    return 1;
}


//create a null matrix of size size*size
double *generate_matrixNull(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = 0;
    }

    return matrix;
}

// create a matrix column of size size
double *generate_matrixColumn(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size);
    srand(1);

    for (i = 0; i < size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}


// print a matrix column named matrix of size size
void print_matrixColumn(const char *name, double *matrix, int size)
{
    int i, j;

    for (i = 0; i < size; i++)
    {
        printf("%f ", matrix[i]);
        printf("\n");
    }
}

// calculate the norm 2 of a matrix named matrix of size n*n
float norm2(double *matrix, int n){
	float norm=0;
	for (int j=0;j<n;j++){
		norm = norm + pow(matrix[j],2);
	}
	norm = sqrt(norm);
	return norm;
}

// calculate the vectorial product between the vector v1 and v2 of size n
double vectorialProduct(double *v1, double *v2, int n){
	double res = 0;
	for (int i=0;i<n;i++){
		res = res + v1[i]*v2[i];
	}
	return res;
}


// get the column indice of a matrix named matrix of size n*n
double *getColumn(double *matrix, int indice, int n){
	double *res=generate_matrixColumn(n);
	for (int i=0;i<n;i++){
		res[i]=0;
	}
	for (int i=0;i<n;i++) {
		res[i] = matrix[indice + i*n];
	}
	return res;
}

// multiply a column of a matrix of size n by a double named R
double *multiplicationScalar(double R, double *colQ, int n){
	double *res = generate_matrixColumn(n);
	for (int i=0;i<n;i++){
		res[i] = R*colQ[i];
	}
	return res;
}


// calculate the result of a vector (tmp of size n) division by a double named norm
double *calcDiv(double *tmp, double norm, int n){
	double *res = generate_matrixColumn(n);
	for (int i=0;i<n;i++){
		res[i] = tmp[i]/norm;
	}
	return res;
}

// calculate the result of a soustraction between two matrix colums of size n
double *soustractionMatrix(double *colp, double *firststep, int n){
	double *res = generate_matrixColumn(n);
	for (int i=0;i<n;i++){
		res[i] = colp[i] - firststep[i];
	}
	return res;
}

// transpose a matrix named matrix of size n*n
double *transpose(double *matrix, int n){
	double *transposer = generate_matrix(n);
	for (int i=0;i<n*n;i++){
		transposer[i] = matrix[i*n%(n*n-1)];
	}
	transposer[n*n-1] = matrix[n*n-1];
	return transposer;
}

// calculate the matricial product between a matrix a and a matrix b of size n*n
double *matricialProduct(double *a,double *b, int n){
	double *res = generate_matrixNull(n);
	for (int I=0; I<n; I++) {
    	for (int J=0; J<n; J++){
        	res[I*n+J]=0;
          	for (int K=0; K<n; K++) {
               res[I*n+J] += a[I*n+K]*b[K*n + J];
          	}
		}
	}
	return res;
}

// calculate the QR decomposition of a matrix A of size n*n and solve the system QRX = B
double *QRCalc(int n, double *A, double *B){
	// initialise variables
	double *Q = generate_matrixNull(n);
	double *R = generate_matrix(n);
	double *p = generate_matrixNull(n);
	double *restest = generate_matrixNull(n);
	double norm = 0;
	double *tmp = getColumn(A,0,n);
	norm = norm2(tmp,n);
	R[0] = norm;
	restest[0] = norm;
	int nb_op = 2*n-1;
	// calculate Q and R
	double *col1Q = calcDiv(tmp,norm,n);
	for (int i=0;i<n;i++){
		Q[i + i*(n-1)] = col1Q[i];
	}
	
	nb_op = nb_op + n;
	for (int i=0;i<(n-1);i++){
		double *colA = getColumn(A,i+1,n);
		int acc = i + 1;
		for (int j=0;j<n;j++){
			p[j + acc] = colA[j];
			acc = acc + (n - 1);
		}
		for (int k=0;k<(i+1);k++){
			double *colQ = getColumn(Q,k,n);
			R[k + i + 1] = vectorialProduct(colA,colQ,n);
			if (k==0){
				restest[i+1] = R[k + i + 1];
			} else {
				restest[k*n + i + 1] = R[k + i +1];
			}
			nb_op = nb_op + 2*n-1;
			double *premiereetape = multiplicationScalar(R[k + i + 1],colQ,n);
			double *colp = getColumn(p,i+1,n);
			double *deuxiemeetape = soustractionMatrix(colp,premiereetape,n);
			for (int m=0;m<n;m++){
				colp[m] = deuxiemeetape[m];
			}
			acc = i+1;
			for (int j=0;j<n;j++){
				p[j + acc] = colp[j];
				acc = acc + (n - 1);
			}
			nb_op = nb_op + 3*n;
		}
		double *colp = getColumn(p,i+1,n);
		R[i*2+3] = norm2(colp,n);
		restest[i*(n+1)+(n+1)] = norm2(colp,n);
    	nb_op = nb_op + 2*n-1;
    	double *colq = getColumn(Q,i+1,n);
    	double *premierres = calcDiv(colp,R[i*2 + 3],n);
    	for (int i=0;i<n;i++){
    		colq[i]=premierres[i];
    	}
    	int acc2 = i+1;
    	for (int i=0;i<n;i++){
    		Q[i + acc2] = colq[i];
    		acc2 = acc2 + (n - 1);
    	}
    	nb_op = nb_op + n;
    }
	//calculate the value of X
	double *X = generate_matrixNull(n);
	double *Qtransposer = transpose(Q,n);
	double *QtransposerB = matricialProduct(Qtransposer,B,n);
	for (int j=n-1;j>=0;--j){
		for (int i=n-1;i>=0;--i){
			double sum = 0;
			for (int k=n-1;k>=0;--k){
				sum = sum + restest[k + i*n]*X[k*n + j];
			}
			if (restest[j + i*n] != 0){
				X[j + i*n] = (QtransposerB[j + i*n] - sum)/restest[j + i*n];
			}
		}
	}
	return X;
}


    void main(int argc, char *argv[])
    {

        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size);
        aref = generate_matrix(size);        
        b = generate_matrix(size);
        bref = generate_matrix(size);

        // Using MKL to solve the system
        MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        tStart = clock();    
        double *X = QRCalc(atoi(argv[1]),a,b);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        if (check_result(bref,X,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        
        print_matrix("X", X, size);
        print_matrix("Xref", bref, size);
    }
