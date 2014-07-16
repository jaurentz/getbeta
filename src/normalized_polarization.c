//** this function will take a vector function. Normalize it, and finally return the polarization of it.*****

//this function call relies upon the **trapz** function and also the **polarization** function. 

#include <stdlib.h>
#include <math.h>
#include <getbeta.h>

void normalize_polarize(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral){
	//cast ii as int to use in for loop. 
	int ii;
 
	//defining the vector psi_hat and dynamically allocating memory for it. 
	double * psi_hat;
	psi_hat = (double*)malloc(num_points*sizeof(double));

	//multiplying the two functions in order to get one function(vector). 
	for(ii=0;ii<num_points;ii++){
		psi_hat[ii] = function_1[ii] * function_1[ii];
	}

	//taking the integral. 
	trapz(psi_hat, left_endpoint, right_endpoint, num_points, integral);

	//finding the norm of that integral. 
	double norm = pow(*integral,0.5); 

	//normalizing the function in question. 
	for(ii=0;ii<num_points;ii++){
		psi_hat[ii] = function_1[ii] / norm;
	}

	//call to polarization function to take normalized vector and return the polarization. e.g. <psi_hat|psi_hat*x>
	polarization(&psi_hat[0], left_endpoint, right_endpoint, num_points, integral);


}



