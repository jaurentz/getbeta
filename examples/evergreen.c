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


int main(){

	struct HamData *HD;
	HD = (struct HamData*)malloc(sizeof(struct HamData));
	HD-> num_points = 100.0; //num_points
	HD-> right_endpoint = 10.0; //right_endpoint formerly known as 'a'
	HD-> left_endpoint = -10.0; //left_endpoint formerly known as 'b'


	double integral = 0; 
	int num_points = HD-> num_points;
	int ii; 
	int jj;
	int kk; 
	int ll; 


	double e_field_max = 2.0; //upper & lower domain of the electric field. 
	double e_field_min = -2.0;

	int row = 0; //setting the row value for the chebyshev differentiation matrix. 
	int col = 0; //setting the column value for the chebyshev differentiation matrix. 

	double sum; 
	double green;
	double BETA = 0;

	//defining the constant alpha and num_chebpoints
	int num_chebpoints = 11; // this should be an odd number so there is a point at zero. 
	double alpha = (pow(HD -> num_points + 1.0, 2.0)) / (pow(HD -> right_endpoint - HD -> left_endpoint,2.0)); 
	double omega = 1.0;
	int flag = 1.0; 
	int tag = 0;


//************************************************************************


	double *main_diag; //main diagonal for the lapack solver
	main_diag = (double*)malloc(num_points*sizeof(double));

	double *sub_diag; //sub diag for the lapack solver(since this matrix is symetric the sub and super diagonals are the same. 
	sub_diag = (double*)malloc(num_points*sizeof(double));

	double *Z; // identity matrix for the lapacke solver 
	Z = (double*)malloc(num_points*num_points*sizeof(double));

	double *psi; //a place to store the first eigenvector.  
	psi = (double*)malloc(num_points*sizeof(double));

	double *cheb_points;//place to store the chebpoints 
	cheb_points = (double*)malloc(num_chebpoints*sizeof(double)); 

	double *elec_field; //place to store the electric field points which correspond to the chebpoints. 
	elec_field = (double*)malloc(num_chebpoints*sizeof(double)); 

	double *diff; //place to store the value of 
	diff = (double*)malloc(sizeof(double));

	double *polarization; //the polarization vector. 
	polarization = (double*)malloc(num_chebpoints*sizeof(double)); 

	//this assigns the electric field array chebpoints on the user defined e_min & e_max so that was it can vary the electric field using the chebyshev points. 
	cheb_array(num_chebpoints, &elec_field[0], e_field_max, e_field_min); 


//************************************************************************


	//creating giant for loop in order to vary the electric field each time. **************
	for(kk=0;kk<num_chebpoints;kk++){

	//*****creating the hamiltonian*****
 
	lap_explicit(HD -> num_points, &main_diag[0], &sub_diag[0], alpha, 0);

	x_squared(HD -> num_points, HD-> left_endpoint, HD-> right_endpoint, &main_diag[0], omega, 1);

	E_field( HD -> num_points, HD -> left_endpoint, HD-> right_endpoint, &main_diag[0], elec_field[kk] , 1);

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



	//Call 1-800-LA-PACKE for a free consulatation of eigenvalues and eigenvectors. 

	LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I', num_points, main_diag, sub_diag, Z, num_points); 


	//this function will take the first eigenvector from the lapacke solver and compute the polarization using the various functions that I've defined. 
	normalize_polarize(&Z[0] ,HD -> left_endpoint ,HD-> right_endpoint ,num_points ,&integral);

	//this will store the 'integral' value (e.g. the polarization in the polarization vector since the normalize_polarize function is a void pointer. So, this will add the data to a vector which can be manipulated and called later on. 
	polarization[kk] = integral; 


}
//end of the giant for loop creating the polarization vector.  

//************************************************************************

	//this takes the polarization vector and will take the inner product of it with the spectral differentiation matrix at zero(provided the number of chebpoints is odd.) This will return the value of BETA. 
	BETA = beta(num_chebpoints, &elec_field[0], &polarization[0]);

	printf("This is beta %E \n", BETA); 


	//freeing the memory of the vectors. 
	free(main_diag);
	free(sub_diag);
	free(cheb_points); 
	free(diff); 

	return 0; 
} 





















