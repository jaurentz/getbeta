//this function will add to the main diagonal the electric field operator which is varied using the alpha. The alpha would be a vector full of chebpoints in order to get different polarization values. 

//the flag must not be zero(0) in order for this operator to add to the main diagonal. 


#include <stdlib.h>
#include <math.h>


void E_field(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag){

	double interval_length = (right_endpoint - left_endpoint)/( num_points +1.0);
	int ii; 

	if(flag == 0){

	for(ii=0;ii<num_points;ii++){
		main_diag[ii] = ( left_endpoint + (ii+1.0)*interval_length) * alpha;		
	}

	}	

	else{ 
	for(ii=0;ii<num_points;ii++){
		main_diag[ii] += ( left_endpoint + (ii+1.0)*interval_length) * alpha;		
	}

	}	

	

}
