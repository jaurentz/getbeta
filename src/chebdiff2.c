//this function will create the second derivative spectral differentiation matrix. 

#include <stdlib.h>
#include <math.h>
#include <getbeta.h>


double chebdiff2(int row, int col, int num_chebpoints, double *cheb_points){


	int middle = num_chebpoints / 2.0;
	int ii;
	int jj = col; 
	double sum = 0; 

	//for(jj=0;jj<num_chebpoints;jj++){
		for(ii=0;ii<num_chebpoints;ii++){

		sum += chebdiff(middle, ii, num_chebpoints, cheb_points) * chebdiff(ii, jj, num_chebpoints, cheb_points); 
	
		}
	
	return sum; 

}

