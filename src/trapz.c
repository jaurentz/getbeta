//this function takes an array and will integrate using a trapezoidal method. It assumes the endpoints are zero. 
#include<getbeta.h>

void trapz(double * function_values, double left_endpoint, double right_endpoint, int num_points, double *integral){   
	
	int ii;	
			
	double constant = (right_endpoint - left_endpoint)/(2.0*(num_points +1.0)) ;

	for(ii=0;ii<num_points;ii++){
		*integral += constant * 2.0 * (function_values[ii]);
	} 
}


