#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <getbeta.h> 

struct HamData
{
	int n;
	double a;
	double b;
};
 
void HamOp(double*X,double*Y,struct HamData*HD);


int main(){

	struct HamData *HD;
	HD = (struct HamData*)malloc(sizeof(struct HamData));
	HD-> n = 3;
	HD-> a = 2.0;
	HD-> b = -2.0;

	int n = HD-> n;

	//allocating dynamic memory
	double *X; //pointer to double array
	X = (double*)malloc(n*n*sizeof(double));
	
	double *Y; //pointer to double array
	Y = (double*)malloc(n*n*sizeof(double));

	double *d; //for the lapack solver
	d = (double*)malloc(n*sizeof(double));

	double *e; //for the lapack solver
	e = (double*)malloc(n*sizeof(double));

	double *Z; //for the lapack solver 
	Z = (double*)malloc(n*n*sizeof(double));

	// for loop to change X  
	int jj;	
	int ii;
	for(ii=0;ii<n;ii++){
		X[ii] = 1.0;
	}
	
	//Assigning the value of Y (zero is a placeholder because it will be computed later.)
	for(ii=0;ii<n;ii++){
		Y[ii] = 0.0;
	}


	//creating the identity matrix which will be used in order to print the Hamiltonian. (H * identity = Hamiltonian) this way we can verify what the Hamiltonian looks like.   
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++){
			if(ii == jj) {
				X[ii*n+jj] = 1.0; 
				//printf("X[%d] = %+1.15e\n",ii,X[ii]);
			}
			else{
				X[ii*n+jj] = 0.0; 
				//printf("X[%d] = %+1.15e\n",jj,X[jj]);
			}
		}
	}
	/*/ printing the identity matrix in one long array 
	printf("\nThis is the new X:\n");
	for(ii=0;ii<n*n;ii++){
		printf("X[%d] = %+1.15e\n",ii,X[ii]);
	}
	printf("\n");			
	
	printf("This is the identity matrix:\n");
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++)
		{ printf("%+1.1e ",X[ii+jj*n]);
		}printf("\n");
	}/*/
		
	//creating the hamiltonian 
	for(ii=0;ii<n;ii++){
		HamOp(&X[ii*n],&Y[ii*n],HD);
	} 
	//printing hamiltonian 
	printf("The is the Hamiltonian:\n");
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++)
		{ printf("%+1.1e ",Y[ii+jj*n]);
		}printf("\n");
	}


	//picking off the main diagonal in order to feed to lapack's eigen-solver 
	printf("Thse are the diagonal elements\n"); 
	for(ii=0;ii<n;ii++){
		d[ii] = Y[ii*(n+1)];
	}
	/*/printing the main diagonal elements. 
	for(ii=0;ii<n;ii++){
	printf("%+1.1e",d[ii]); 
	}printf("\n");/*/

	//picking the off-diagonal elements in order to feed to lapack's eigen-solver
	printf("These are the off-main diagonal elements");
	for(ii=0;ii<n-1;ii++){
		e[ii] = Y[ii*(n+1)+1];
	}printf("\n");

	/*/printing the off-diagonal elements. 
	for(ii=0;ii<n-1;ii++){
	printf("%+1.1e",e[ii]); 
	}printf("\n");/*/

	//identity matrix(Z) Needed for lapack's eigen-solver(otherwise you won't get the correct values)    
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++){
			if(ii == jj) {
				Z[ii*n+jj] = 1.0; 
				//printf("X[%d] = %+1.15e\n",ii,X[ii]);
			}
			else{
				Z[ii*n+jj] = 0.0; 
				//printf("X[%d] = %+1.15e\n",jj,X[jj]);
			}
		}
	}
	/*/printing the identity matrix. 
	printf("This is the new identity matrix:\n");
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++)
		{ printf("%+1.1e ",Z[ii+jj*n]);
		}printf("\n");
	}/*/

	//Call 1-800-LA-PACKE for a free consulatation of eigenvalues and eigenvectors. 

	LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I', n, d ,e ,Z, n); 



	
	/*
	printf("\nThis is the new Y:\n");
	for(ii=0;ii<n;ii++){
		printf("Y[%d] = %+1.15e\n",ii,Y[ii]);
	}
	printf("\n");/*/
	
	
	//freeing memory 
	free(X);
	free(Y);
	free(d);
	free(e);
	free(Z);
	


	 
    return 0;
}


void HamOp(double*X,double*Y,struct HamData*HD){


	double h = (HD-> a- HD-> b)/(HD->n+1.0);
	
	int flag =  0;
	
	//constant for the laplacian e.g. the 1/h^2 term 
	double alpha = (1.0/pow(h,2.0));

	lapmult(alpha ,HD-> n , X, Y, flag);
	
	//constant for the potential 
	double beta = 1.0;

	flag  = 1;

	potential(beta ,HD->n ,HD-> a ,HD-> b , X, Y, flag);

 
}









