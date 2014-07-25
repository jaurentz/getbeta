//this function will give back the proper of element of the spectral differentiation matrix defined by the row and col. This is taken from 'Spectral Methods in MATLAB' by Trefethen p.53. 

//this function needs to be pass cheb_points. 


#include <stdlib.h>
#include <math.h>
#include <getbeta.h>


double chebdiff(int row, int col, int num_chebpoints, double *cheb_points){

	int N = num_chebpoints-1; 
	//int row = 3; 
	//int col = 3;
	double C_i;
	double C_j ;
	double C ;

	double diff;

	if (row == 0 && col == 0){
		diff = 1.0/6.0*(2.0*pow(N,2.0)+1.0);
	}

	else if(row == N && col == N){
		diff = -(1.0/6.0*(2.0*pow(N,2.0)+1.0));
	}
	
	else if(row == col){
		diff = -(cheb_points[row])/ (2.0 * (1.0 - pow(cheb_points[row],2.0)));
	}
	
	else{	
		if(row == 0 || row == N){C_i = 2.0;} // the OR statement in C is "||"!! 
		else{C_i = 1.0;}
		if(col == 0 || col == N){C_j = 2.0;}
		else{C_j = 1.0;}

		C = C_i / C_j; 

		diff = C * (pow(-1.0, row + col) / (cheb_points[row] - cheb_points[col]));
	}

	return diff;

}
