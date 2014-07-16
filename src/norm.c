//this function will compute the NORM of a vector. You would have to divide the original vector (element-wise) by the value this function returns in order to NORMALIZE the vector. 

//This function relies on calling the **TRAPZ** function. 
#include<stdlib.h>
#include<math.h>
#include<getbeta.h>


void norm(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral){
	//cast ii as int so use in for loop. 
	int ii; 
	//creating the pointer psi_hat and dynamically allocating memory for it. 
	double * psi_hat;
	psi_hat = (double*)malloc(num_points*sizeof(double));
	//multiplying the two vectors(functions). 
	for(ii=0;ii<num_points;ii++){
		psi_hat[ii] = function_1[ii] * function_1[ii];
	}
	//calling trapz function to integrate the new psi_hat function. 
	trapz(psi_hat, left_endpoint, right_endpoint, num_points, integral);

	//computes the norm by taking the square root of the computed integral. 
	double norm = pow(*integral,0.5); 

	/*this section would normalize the vector, if you wanted to do something like that. 
	for(ii=0;ii<num_points;ii++){
		psi_hat[ii] = function_1[ii] / norm;
	}/*/

	free(psi_hat);

}

