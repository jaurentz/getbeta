#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <getbeta.h> 

void lap_explicit(int num_points, double *main_diag, double *sub_diag, double alpha, int flag);

void x_squared(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag);

void E_field(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag);

//void normalize_polarize(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral);

//void polarization(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral);

struct HamData
{
	int num_points; //num_points
	double right_endpoint; //right_endpoint
	double left_endpoint; //left_endpoint 
};
 
void HamOp(double*X,double*Y,struct HamData*HD);


int main(){

	struct HamData *HD;
	HD = (struct HamData*)malloc(sizeof(struct HamData));
	HD-> num_points = 10; //num_points
	HD-> right_endpoint = 2.0; //right_endpoint formerly known as 'a'
	HD-> left_endpoint = -2.0; //left_endpoint formerly known as 'b'

	double integral = 0; 
	int num_points = HD-> num_points;
	int ii; 
	int jj;
	int kk; 
	int ll; 
	double  elec_field[3] ={-1,0,1};
 
	double p[3] = {0,0,0};
	int m = 3;  

	//defining the constant alpha 

	double alpha = 1.0; 
	int flag = 1.0; 

	double *main_diag; //for the lapack solver
	main_diag = (double*)malloc(num_points*sizeof(double));

	double *sub_diag; //for the lapack solver
	sub_diag = (double*)malloc(num_points*sizeof(double));

	double *Z; // identity matrix for the lapacke solver 
	Z = (double*)malloc(num_points*num_points*sizeof(double));

	double *psi; //a place to store the first eigenvector.  
	psi = (double*)malloc(num_points*sizeof(double));




	for(kk=0;kk<m;kk++){

	//*****creating the hamiltonian*****
 
	lap_explicit(HD -> num_points, &main_diag[0], &sub_diag[0], alpha, flag);

	x_squared(HD -> num_points, HD-> left_endpoint, HD-> right_endpoint, &main_diag[0], alpha, flag);
	
	//for(jj=0;jj<)
	E_field( HD -> num_points, HD -> left_endpoint, HD-> right_endpoint, &main_diag[0], elec_field[kk] , flag);



	//identity matrix(Z) Needed for lapack's eigen-solver(otherwise you won't get the correct values)    
	for(ii=0;ii<num_points;ii++)
	{
		for(jj=0;jj<num_points;jj++){
			if(ii == jj) {
				Z[ii*num_points+jj] = 1.0; 
				//printf("X[%d] = %+1.15e\n",ii,X[ii]);
			}
			else{
				Z[ii*num_points+jj] = 0.0; 
				//printf("X[%d] = %+1.15e\n",jj,X[jj]);
			}
		}
	}



	//Call 1-800-LA-PACKE for a free consulatation of eigenvalues and eigenvectors. 

	LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I', num_points, main_diag ,sub_diag ,Z, num_points); 

	/*
	pick off the first eigenvector eg the ground state of the Hamiltonian. 
	for(ii=0;ii<num_points;ii++){
		//psi[ii] = Z[ii]; 




	printf("\nThis is Psi:\n");
	for(ii=0;ii<num_points;ii++){
		printf("Psi[%d] = %+1.15e\n",ii,psi[ii]);
	}
	printf("\n");/*/




	normalize_polarize(&Z[0] ,HD -> left_endpoint ,HD-> right_endpoint ,num_points ,&integral);

	printf("The polarization: %1.1e \n", integral);	


}


	free(main_diag);
	free(sub_diag);

	return 0; 
} 

void lap_explicit(int num_points, double *main_diag, double *sub_diag, double alpha, int flag){

	int ii; 

	if(flag == 0){
	for(ii=0;ii<num_points;ii++){
		main_diag[ii]= 2.0;
	}

	for(ii=0;ii<num_points - 1;ii++){
		sub_diag[ii]= -1.0;
	}

	}

	else{
	for(ii=0;ii<num_points;ii++){
		main_diag[ii] += 2.0;
	}

	for(ii=0;ii<num_points - 1;ii++){
		sub_diag[ii] = -1.0;
	}

	}


}




void x_squared(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag){
	
	double interval_length = ( right_endpoint -  left_endpoint)/( num_points +1.0);
	int ii; 	

	if(flag == 0){
	for(ii=0;ii<num_points;ii++){
		main_diag[ii] = (pow(( left_endpoint + (ii+1.0)*interval_length),2)) * alpha;		
	}

	}

	else{
		for(ii=0;ii<num_points;ii++){
		main_diag[ii] += (pow(( left_endpoint + (ii+1.0)*interval_length),2)) * alpha;		
	}

	}

}


void E_field(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag){

	double interval_length = (right_endpoint - left_endpoint)/( num_points +1.0);
	int ii; 

	if(flag == 0){

	for(ii=0;ii<num_points;ii++){
		main_diag[ii] = ( left_endpoint + (ii+1.0)*interval_length) * alpha;		
	}

	}	

	else{ 
	for(ii=0;ii<num_points;ii++){
		main_diag[ii] += ( left_endpoint + (ii+1.0)*interval_length) * alpha;		
	}

	}	

	

}


