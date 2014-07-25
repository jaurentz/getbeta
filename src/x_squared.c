//this function will add the x-squared operator to the main-diagonal array. 

//The flag must be 1(or anything but zero(0) ) to add to the main 

//alpha here could be an arbitrary scalar. 

#include <stdlib.h>
#include <math.h>

void x_squared(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag){
	
	double interval_length = ( right_endpoint -  left_endpoint) / ( num_points + 1.0);
	int ii; 	

	if(flag == 0){
	for(ii=0;ii<num_points;ii++){
		main_diag[ii] = (pow(( left_endpoint + (ii+1.0)*interval_length),2.0)) * alpha;		
	}

	}

	else{
		for(ii=0;ii<num_points;ii++){
		main_diag[ii] += (pow(( left_endpoint + (ii+1.0)*interval_length),2.0)) * alpha;		
	}

	}

}

