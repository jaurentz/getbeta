//this function computes the polarization. This means it will take a psi^2*x and integrate over the the interval. 

//function need **x_operator** & **inner_product**  

#include <stdlib.h>
#include<getbeta.h>

void polarization(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral){
		
	int flag = 1;
	int ii; 
	double beta = 1.0; 
	//defining the size of the interval(s). 
	double h = (right_endpoint - left_endpoint)/(num_points+1.0); 

	//allocating memory for new vector for X * vector 
	double * psi_hat; 
	psi_hat = (double*)malloc(num_points*sizeof(double));

	X_operator(beta, num_points, right_endpoint, left_endpoint, psi_hat, function_1, flag); 

	inner_product(&function_1[0], &psi_hat[0], left_endpoint, right_endpoint, num_points, integral);   

	free(psi_hat);

}

