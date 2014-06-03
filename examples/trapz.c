#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void trapz(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double *integral);   

void inner_product(double * function_1, double * function_2, double left_endpoint, double right_endpoint, int interior_points, double *integral);

void polarization(double * function_1,double left_endpoint, double right_endpoint, int interior_points, double * integral); 

void X_operator(double beta, int n, double a, double b, double *psi_hat, double *function_1, int flag);

int main(){

	double integral = 0.0; 

	int n = 8;
	double Z[n];


	double b = 3.0; 
	double a = 0.0;

	int ii; 
	
	for(ii=0;ii<n;ii++){
		Z[ii] = 1.0;
	}cd

	/*for(ii=0;ii<n;ii++){
	printf("Z: %+1.1e",Z[ii]); 
	}printf("\n");/*/

	//trapz(&Z[0],a,b,n,&integral);
	
	//inner_product(&Z[0],&Z[0],a,b,n,&integral);

	//printf("The inner product: %1.1e \n", integral);

	polarization(&Z[0], a, b, n, &integral);

	printf("The polarization: %1.1e \n", integral);	

	//printf("Integral: %+1.1e \n",integral);
	
	//free(Z);
	



	return 0;
}



void trapz(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double *integral){   
	
	int ii;	
			
	double constant = (right_endpoint - left_endpoint)/(2.0*(interior_points +1.0)) ;

	for(ii=0;ii<interior_points;ii++){
		*integral += constant * 2.0 * (function_values[ii]);
	} 
}

void inner_product(double * function_1, double * function_2, double left_endpoint, double right_endpoint, int interior_points, double * integral){
	
	int ii; 

	double * psi_hat;
	psi_hat = (double*)malloc(interior_points*sizeof(double));

	for(ii=0;ii<interior_points;ii++){
		psi_hat[ii] = function_1[ii] * function_2[ii];
	}

	trapz(psi_hat, left_endpoint, right_endpoint, interior_points, integral);

	free(psi_hat);

}

void polarization(double * function_1, double left_endpoint, double right_endpoint, int interior_points, double * integral){
		
	int flag = 1;
	int ii; 
	double beta = 1.0; 
	//defining the size of the interval(s). 
	double h = (right_endpoint - left_endpoint)/(interior_points+1.0); 

	//allocating memory for new vector for X * vector 
	double * psi_hat; 
	psi_hat = (double*)malloc(interior_points*sizeof(double));

	X_operator(beta, interior_points, right_endpoint, left_endpoint, psi_hat, function_1, flag); 

	inner_product(&function_1[0], &psi_hat[0], left_endpoint, right_endpoint, interior_points, integral);   

	free(psi_hat);

}


void X_operator(double beta, int n, double a, double b, double* psi_hat, double* function_1, int flag){

	int jj;
	double h = (a - b)/(n+1.0);

	if(flag == 0){	
	//jj+1 because it picks the first interval + 1
	for(jj=0;jj<n;jj++){
		psi_hat[jj] = ( b + (jj+1.0)*h) * (function_1[jj])*beta; 
	}
	
	}
	else{  //if else do this
	for(jj=0;jj<n;jj++){
		psi_hat[jj] += ( b +(jj+1.0)*h) * (function_1[jj])*beta;
	}

	}
 

}




