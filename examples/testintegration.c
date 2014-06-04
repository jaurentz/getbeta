#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getbeta.h> 


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

	trapz(&Z[0],a,b,n,&integral);
	
	printf("The inner product: &1.1e \n", integral);	

	//inner_product(&Z[0],&Z[0],a,b,n,&integral);

	//printf("The inner product: %1.1e \n", integral);

	//polarization(&Z[0], a, b, n, &integral);

	//printf("The polarization: %1.1e \n", integral);	

	//printf("Integral: %+1.1e \n",integral);
	
	//free(Z);
	



	return 0;
}







