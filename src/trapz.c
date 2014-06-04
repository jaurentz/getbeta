//this function takes an array and will integrate using a trapezoidal method. It assumes the endpoints are zero. 


void trapz(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double *integral){   
	
	int ii;	
			
	double constant = (right_endpoint - left_endpoint)/(2.0*(interior_points +1.0)) ;

	for(ii=0;ii<interior_points;ii++){
		*integral += constant * 2.0 * (function_values[ii]);
	} 
}


