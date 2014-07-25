//this function should give the chebpoints for an arbitrary interval [e_field min, e_field_max] and will store them in the vector cheb_points[ii]


#include <stdlib.h>
#include <math.h>


void cheb_array(int num_chebpoints, double *cheb_points, double e_field_max, double e_field_min){

	int ii;

	for(ii=0;ii<num_chebpoints;ii++){

	cheb_points[ii] = -1.0*(((e_field_max - e_field_min))/(2.0)) * sin(M_PI*((ii/(num_chebpoints-1.0)) - 0.5)) + ((e_field_min + e_field_max)/(2.0));	
	}

}

