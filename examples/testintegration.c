#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getbeta.h> 


int main(){

	double integral = 0.0; 

	int num_points = 10000;
	double Z[num_points];
	double Q[num_points];


	double right_endpoint = 3.0; 
	double left_endpoint = 0.0;

	int ii; 
	
	//for(ii=0;ii<n;ii++){
	//	Z[ii] = 1.0;
	//}

	int jj;
	double interval_length = (right_endpoint - left_endpoint)/(num_points + 1.0);

	//printf("This is h %1.1e \n", h);
	
	//creating a sine function 
	//ii+1 because it picks the first interval + 1
	for(ii=0;ii<num_points;ii++){
		Z[ii] = sin((2.0 * M_PI)/(right_endpoint - left_endpoint)* (( left_endpoint + (ii+1.0)* interval_length)- left_endpoint) ) ; 
	}

	/*
	for(ii=0;ii<n;ii++){
	printf("Z: %+1.1e",Z[ii]); 
	}printf("\n");/*/ 

	//creating a different function that is orthogonal to the first sine function.
	//ii+1 because it picks the first interval + 1
	for(ii=0;ii<num_points;ii++){
		Q[ii] = sin((4.0 * M_PI)/(right_endpoint - left_endpoint)* (( left_endpoint + (ii+1.0)* interval_length)- left_endpoint) ) ; 
	}



	//trapz(&Z[0], left_endpoint, right_endpoint, num_points, &integral);
		
	//printf("The integral of sin: %1.1e \n", integral);	

	//inner_product(&Z[0],&Q[0], left_endpoint, right_endpoint, num_points, &integral);

	//printf("The inner product: %1.1e \n", integral);

	polarization(&Z[0], left_endpoint, right_endpoint, num_points, &integral);

	printf("The polarization: %1.1e \n", integral);	

	//printf("Integral: %+1.1e \n",integral);
	
	//free(Z); you can't unallocate STATIC memory!! 
	



	return 0;
}







