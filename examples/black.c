#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct HamData
{
	int n;
	double a;
	double b;
};
 
void HamOp(double*X,double*Y,struct HamData*HD);



int main(){

	struct HamData *HD;
	HD = (struct HamData*)malloc(sizeof(struct HamData));
	HD-> n = 7.0;
	HD-> a = 2.0;
	HD-> b = -2.0;

	int n = HD-> n;

	//creating dynamic memory
	double *X; //pointer to double array
	X = (double*)malloc(n*sizeof(double));
	
	double *Y; //pointer to double array
	Y = (double*)malloc(n*sizeof(double));


	
	//computing the size of the interval
	double h = (HD-> a-HD-> b)/(n+1.0);

	// for loop to change X
	int ii;
	for(ii=0;ii<n;ii++){
		X[ii] = 2.0;
	}
	
	//for loop to change Y
	for(ii=0;ii<n;ii++){
		Y[ii] = 0.0;
	}

	HamOp(X,Y,HD); 
	
	printf("\nThis is the new Y:\n");
	for(ii=0;ii<n;ii++){
		printf("Y[%d] = %+1.15e\n",ii,Y[ii]);
	}
	printf("\n");
	

	//freeing memory 
	free(X);
	free(Y);
	


	 
    return 0;
}


void HamOp(double*X,double*Y,struct HamData*HD){

	double h = (HD-> a- HD-> b)/(HD->n+1.0);

	int ii;
	//computing Y
	for(ii=1;ii<HD-> n-1;ii++){
		Y[ii] = -X[ii-1]+2*X[ii]-X[ii+1];
	}
	//top row of laplacian 
	Y[0] = 2*X[0]-X[1];
	
	//bottom row of laplacian 
	Y[HD->n-1] = 2*X[HD->n-1]-X[HD->n-2];

	
	//multiplying all terms by 1/h^2 & then adding the potential operator elements
	for(ii=0;ii<HD->n;ii++){
		Y[ii] = Y[ii]*(1.0/pow(h,2.0)) + pow((HD->b + (ii)*h),2) * X[ii]; 
	}
	
	/*/creating the potential operator 
	int ii; 
	for(ii=0;ii<n;ii++){
		Y[ii] = pow((b + (ii)*h),2) * X[ii]; 
	}/*/

 
}





