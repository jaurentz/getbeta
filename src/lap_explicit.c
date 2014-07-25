//this function will create the main diagonal of the laplacian operator for the lapacke eigen solver. It will also assign the super and sub diagonal entries and place them in an array. 

//It will also multiply by the laplacian constant alpha which is computed outside of the function call. 

//the flag should be 0 because you want to overwrite what is on the main diag and sub/super diags each time rather than adding to them(a running sum). 

#include <stdlib.h>
#include <math.h>


void lap_explicit(int num_points, double *main_diag, double *sub_diag, double alpha, int flag){

	int ii; 
	//writing the value of 2 to main diag. and -1 to sub/super-diag. 
	if(flag == 0){
	for(ii=0;ii<num_points;ii++){
		main_diag[ii]= 2.0 * alpha;
	}

	for(ii=0;ii<num_points - 1;ii++){
		sub_diag[ii]= -1.0 * alpha;
	}

	}

	else{
	for(ii=0;ii<num_points;ii++){
		main_diag[ii] += 2.0 * alpha;
	}

	for(ii=0;ii<num_points - 1;ii++){
		sub_diag[ii] = -1.0 * alpha;
	}

	}

}
