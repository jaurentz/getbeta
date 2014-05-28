#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void lapmult(int n,double*X,double*Y);

int main(){

	
	// allocate dynamic memory
	int n = 10;
	double *X; //pointer to double array
	X = (double*)malloc(n*sizeof(double));

	double *Y; //pointer to double array
	Y = (double*)malloc(n*sizeof(double));

	// for loop to change X
	int ii;
	for(ii=0;ii<n;ii++){
		X[ii] = 2;
	}
	//for loop to change Y
	for(ii=0;ii<n;ii++){
		Y[ii] = 0;
	}

	lapmult(n,X,Y);

	// for loop to print Y
	
	printf("\nThis is Y:\n");
	for(ii=0;ii<n;ii++){
		printf("Y[%d] = %+1.15e\n",ii,Y[ii]);
	}
	printf("\n");

	// free memory
	free(X);
	free(Y);

	// return value
	return 0;

}


void lapmult(int n,double*X,double*Y){

	int ii;
	//computing Y
	for(ii=1;ii<n-1;ii++){
		Y[ii] = -X[ii-1]+2*X[ii]-X[ii+1];
	}
	//top row of laplacian 
	Y[0] = 2*X[0]-X[1];
	
	//bottom row of laplacian 
	Y[n-1] = 2*X[n-1]-X[n-2];
}
