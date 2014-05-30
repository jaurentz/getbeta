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

void lapmult(double alpha, int n, double*X, double*Y, int flag);	

void potential(double beta, int n, double a, double b, double*X, double*Y, int flag); 


int main(){

	struct HamData *HD;
	HD = (struct HamData*)malloc(sizeof(struct HamData));
	HD-> n = 7;
	HD-> a = 2.0;
	HD-> b = -2.0;

	int n = HD-> n;

	//creating dynamic memory
	double *X; //pointer to double array
	X = (double*)malloc(n*n*sizeof(double));
	
	double *Y; //pointer to double array
	Y = (double*)malloc(n*n*sizeof(double));

	// for loop to change X
	int jj;	
	int ii;
	for(ii=0;ii<n;ii++){
		X[ii] = 1.0;
	}
	
	//for loop to change Y
	for(ii=0;ii<n;ii++){
		Y[ii] = 0.0;
	}


	//creating the identity matrix     
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++){
			if(ii == jj) {
				X[ii*n+jj] = 1.0; 
				//printf("X[%d] = %+1.15e\n",ii,X[ii]);
			}
			else{
				X[ii*n+jj] = 0.0; 
				//printf("X[%d] = %+1.15e\n",jj,X[jj]);
			}
		}
	}
	/*/ printing the identity matrix in one long array 
	printf("\nThis is the new X:\n");
	for(ii=0;ii<n*n;ii++){
		printf("X[%d] = %+1.15e\n",ii,X[ii]);
	}
	printf("\n");			
	/*/
	printf("This is the identity matrix:\n");
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++)
		{ printf("%+1.1e ",X[ii+jj*n]);
		}printf("\n");
	}
		
	//for loop to print that ham
	for(ii=0;ii<n;ii++){
		HamOp(&X[ii*n],&Y[ii*n],HD);
	} 

	printf("Home is where the cured pig is(e.g. HAM):\n");
	for(ii=0;ii<n;ii++)
	{
		for(jj=0;jj<n;jj++)
		{ printf("%+1.1e ",Y[ii+jj*n]);
		}printf("\n");
	}
		
	/*
	printf("\nThis is the new Y:\n");
	for(ii=0;ii<n;ii++){
		printf("Y[%d] = %+1.15e\n",ii,Y[ii]);
	}
	printf("\n");/*/
	
	
	//freeing memory 
	free(X);
	free(Y);
	


	 
    return 0;
}


void HamOp(double*X,double*Y,struct HamData*HD){


	double h = (HD-> a- HD-> b)/(HD->n+1.0);
	
	int flag =  1;
	
	//constant for the laplacian e.g. the 1/h^2 term 
	double alpha = (1.0/pow(h,2.0));

	//constant for the potential 
	double beta = 1.0;

	lapmult(alpha ,HD-> n , X, Y, flag);
	
	potential(beta ,HD->n ,HD-> a ,HD-> b , X, Y, flag);

 
}

void lapmult(double alpha ,int n, double*X,double*Y,int flag){

	int ii;

	if(flag == 0){
	
	//computing Y
	for(ii=1;ii< n-1;ii++){
		Y[ii] = (-X[ii-1]+2*X[ii]-X[ii+1])*alpha;
	}
	//top row of lapl
	Y[0] = (2*X[0]-X[1])*alpha;
	
	//bottom row of laplacian 
	Y[n-1] = (2*X[n-1]-X[n-2])*alpha;

	
	//multiplying all terms by 1/h^2
	for(ii=0;ii< n;ii++){
		Y[ii] = Y[ii]*alpha;
	}

	}

	else{//if else do this

	for(ii=1;ii< n-1;ii++){
		Y[ii] += (-X[ii-1]+2*X[ii]-X[ii+1])*alpha;
	}
	
	Y[0] += (2*X[0]-X[1])*alpha;

	Y[n-1] += (2*X[ n-1]-X[n-2])*alpha;

	/* the old way of multiplying by a constant 
	for(ii=0;ii< n ;ii++){
		Y[ii] = Y[ii]*alpha;

	}/*/

	}

	
}

//creating the potential operator 
void potential(double beta, int n, double a, double b, double*X, double*Y,int flag){

	int jj;
	double h = (a - b)/(n+1.0);

	if(flag == 0){	
	
	for(jj=0;jj<n;jj++){
		Y[jj] = (pow(( b + (jj+1.0)*h),2) * X[jj])*beta; 
	}
	
	}
	else{  //if else do this
	for(jj=0;jj<n;jj++){
		Y[jj] += (pow(( b +(jj+1.0)*h),2) * X[jj])*beta;
	}

	}
 

}








