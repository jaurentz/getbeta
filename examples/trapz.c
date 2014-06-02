#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void trapz(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double *integral);   

void trapz_squared(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double *integral);

//void x_and_trapz_squared(

int main(){

	double integral = 0.0; 

	int n = 3;
	double Z[n];


	double b = 2cd.0; 
	double a = -2.0;

	int ii; 
	
	for(ii=0;ii<n;ii++){
		Z[ii] = 10.0;
	}

	for(ii=0;ii<n;ii++){
	printf("Z: %+1.1e",Z[ii]); 
	}printf("\n");

	//trapz(&Z[0],a,b,n,&integral);
	
	trapz_squared(&Z[0],a,b,n,&integral);	

	printf("Integral: %+1.1e \n",integral);
	
	//free(Z);
	//free(psi_hat);



	return 0;
}



void trapz(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double *integral){   
	
	int ii;	
			
	double constant = (right_endpoint - left_endpoint)/(2.0*(interior_points +1.0)) ;

	for(ii=0;ii<interior_points;ii++){
		*integral += constant * 2.0 * (function_values[ii]);
	} 
}

void trapz_squared(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double * integral){
	
	int ii; 

	double * psi_hat;
	psi_hat = (double*)malloc(interior_points*sizeof(double));

	for(ii=0;ii<interior_points;ii++){
		psi_hat[ii] = pow(function_values[ii],2.0);
	}

	trapz(psi_hat, left_endpoint, right_endpoint, interior_points, integral);

	free(psi_hat);

}
