//this function will compute beta by picking off middle row of the chebyshev differentiation matrix and multiplying it by the polarization vector. 

//this is a running sum because this is essentially taking an inner(dot) product. So it is a*a1+b*b1+etc... Hence it is beta ***+=***.

#include <stdlib.h>
#include <math.h>
#include <getbeta.h>


double beta(int num_chebpoints, double *cheb_points, double *integral){

	int ii;
	int row = num_chebpoints/2; 
	double beta = 0.0;
	
	for(ii=0;ii<num_chebpoints;ii++){
		beta += chebdiff2(row, ii, num_chebpoints, cheb_points ) * integral[ii];
	}

	return beta;
}
