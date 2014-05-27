#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){

	// print hello world
	printf("\nhello world!\n\n");

	// declare a variable
	int myage = 13;
	
	// print myage
	printf("\nmyage = %d\n\n",myage);
	
	// declare a double
	double mydouble = (double)5/2;
	
	// print mydouble
	printf("\nmydouble = %1.15e\n\n",mydouble);
	
	// allocate dynamic memory
	int n = 20;
	double *myvec; //pointer to double array
	myvec = (double*)malloc(n*sizeof(double));
	
	// for loop to print myvec
	int ii;
	printf("\nmyvec:\n");
	for(ii=0;ii<n;ii++){
		printf("myvec[%d] = %+1.15e\n",ii,myvec[ii]);
	}
	printf("\n");
	
	// for loop to change myvec
	for(ii=0;ii<n;ii++){
		myvec[ii] = ii*ii;
	}
	
	// for loop to print myvec
	printf("\nmyvec:\n");
	for(ii=0;ii<n;ii+=2){
		printf("myvec[%d] = %+1.15e\n",ii,myvec[ii]);
	}
	printf("\n");
	
	// free memory
	free(myvec);

	// return value
	return 0;

};

