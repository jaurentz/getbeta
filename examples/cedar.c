#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <getbeta.h> 

struct HamData
{
	int num_points; //num_points
	double right_endpoint; //right_endpoint
	double left_endpoint; //left_endpoint 
};
 
void HamOp(double*X,double*Y,struct HamData*HD);


int main(){

	struct HamData *HD;
	HD = (struct HamData*)malloc(sizeof(struct HamData));
	HD-> num_points = 10; //num_points
	HD-> right_endpoint = 2.0; //right_endpoint formerly known as 'a'
	HD-> left_endpoint = -2.0; //left_endpoint formerly known as 'b'

	int num_points = HD-> num_points;

	//allocating dynamic memory
	double *X; //pointer to double array
	X = (double*)malloc(num_points*num_points*sizeof(double));
	
	double *Y; //pointer to double array
	Y = (double*)malloc(num_points*num_points*sizeof(double));

	double *d; //for the lapack solver
	d = (double*)malloc(num_points*sizeof(double));

	double *e; //for the lapack solver
	e = (double*)malloc(num_points*sizeof(double));

	double *Z; //for the lapack solver 
	Z = (double*)malloc(num_points*num_points*sizeof(double));

	// for loop to change X  this could be changed as it will be the input array into the Hamiltonian. 
	int jj;	
	int ii;
	for(ii=0;ii<num_points;ii++){
		X[ii] = 1.0;
	}
	
	//Assigning the value of Y (zero is a placeholder because it will be computed later.)
	for(ii=0;ii<num_points;ii++){
		Y[ii] = 0.0;
	}


	//creating the identity matrix which will be used in order to print the Hamiltonian. (H * identity = Hamiltonian) this way we can verify what the Hamiltonian looks like.   
	for(ii=0;ii<num_points;ii++)
	{
		for(jj=0;jj<num_points;jj++){
			if(ii == jj) {
				X[ii*num_points+jj] = 1.0; //ones down the main diagonal, this is how you pick them off. note: it is ii* num_points+jj 
				//printf("X[%d] = %+1.15e\n",ii,X[ii]);
			}
			else{
				X[ii*num_points + jj] = 0.0; 
				//printf("X[%d] = %+1.15e\n",jj,X[jj]);
			}
		}
	}
	/*/ printing the identity matrix in one long array 
	printf("\nThis is the new X:\n");
	for(ii=0;ii<num_points*num_points;ii++){
		printf("X[%d] = %+1.15e \n",ii,X[ii]);
	}
	printf("\n");			
	
	printf("This is the identity matrix:\n");
	for(ii=0;ii<num_points;ii++)
	{
		for(jj=0;jj<num_points;jj++)
		{ printf("%+1.1e ",X[ii+jj*n]);
		}printf("\n");
	}/*/
		
	//creating the hamiltonian using the function call HamOp
	for(ii=0;ii<num_points;ii++){
		HamOp(&X[ii*num_points],&Y[ii*num_points],HD);
	} 
	//printing hamiltonian 
	printf("The is the Hamiltonian:\n");
	for(ii=0;ii<num_points;ii++)
	{
		for(jj=0;jj<num_points;jj++)
		{ printf("%+1.1e ",Y[ii+jj*num_points]);
		}printf("\n");
	}


	//picking off the main diagonal in order to pass to lapack's eigen-solver 
	printf("These are the diagonal elements\n"); 
	for(ii=0;ii<num_points;ii++){
		d[ii] = Y[ii*(num_points+1)];
	}
	/*/printing the main diagonal elements. 
	for(ii=0;ii<num_points;ii++){
	printf("%+1.1e",d[ii]); 
	}printf("\n");/*/

	//picking the off-diagonal elements in order to pass to lapack's eigen-solver
	printf("These are the off-main diagonal elements");
	for(ii=0;ii<num_points - 1;ii++){
		e[ii] = Y[ii*(num_points +1)+1];
	}printf("\n");

	/*/printing the off-diagonal elements. 
	for(ii=0;ii<n-1;ii++){
	printf("%+1.1e",e[ii]); 
	}printf("\n");/*/

	//identity matrix(Z) Needed for lapack's eigen-solver(otherwise you won't get the correct values)    
	for(ii=0;ii<num_points;ii++)
	{
		for(jj=0;jj<num_points;jj++){
			if(ii == jj) {
				Z[ii*num_points+jj] = 1.0; 
				//printf("X[%d] = %+1.15e\n",ii,X[ii]);
			}
			else{
				Z[ii*num_points+jj] = 0.0; 
				//printf("X[%d] = %+1.15e\n",jj,X[jj]);
			}
		}
	}
	/*/printing the identity matrix. 
	printf("This is the new identity matrix:\n");
	for(ii=0;ii<num_points;ii++)
	{
		for(jj=0;jj<num_points;jj++)
		{ printf("%+1.1e ",Z[ii+jj*num_points]);
		}printf("\n");
	}/*/

	//Call 1-800-LA-PACKE for a free consulatation of eigenvalues and eigenvectors. 

	LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I', num_points, d ,e ,Z, num_points); 

	printf("first eigenvalue. %1.15e \n", d[0]);

	
	/*
	printf("\nThis is the new Y:\n");
	for(ii=0;ii<num_points;ii++){
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


	double interval_length = (HD-> right_endpoint - HD-> left_endpoint)/(HD-> num_points +1.0);
	
	int flag =  0;
	
	//constant for the laplacian e.g. the 1/ interval_length ^2 term 
	double alpha = (1.0/pow(interval_length ,2.0));

	lapmult(alpha ,HD-> num_points , X, Y, flag);
	
	//constant for the potential 
	double beta = 1.0;

	flag  = 1;

	potential(beta ,HD->num_points ,HD-> right_endpoint ,HD-> left_endpoint , X, Y, flag);

 
}









