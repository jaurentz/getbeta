#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <getbeta.h> 

/*
void lap_explicit(int num_points, double *main_diag, double *sub_diag, double alpha, int flag);

void x_squared(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag);

void E_field(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag);

//void normalize_polarize(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral);

//void polarization(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral);

double chebdiff(int row, int col, int num_chebpoints, double *cheb_points);

double chebdiff2(int row, int col, int num_chebpoints, double *chebpoints);

double beta(int row, int col, int num_chebpoints, double *cheb_points, double *polarization);

void cheb_array(int num_chebponts, double *cheb_points, double e_field_max, double e_field_min);
/*/


struct HamData
{
	int num_points; //num_points
	double right_endpoint; //right_endpoint
	double left_endpoint; //left_endpoint 
};

 
//void HamOp(double*X,double*Y,struct HamData*HD);


int main(){

	struct HamData *HD;
	HD = (struct HamData*)malloc(sizeof(struct HamData));
	HD-> num_points = 100.0; //num_points
	HD-> right_endpoint = 10.0; //right_endpoint formerly known as 'a'
	HD-> left_endpoint = -10.0; //left_endpoint formerly known as 'b'

	double integral = 0; 
	int num_points = HD-> num_points;
	int ii; 
	int jj;
	int kk; 
	int ll; 
	//double  elec_field[3] ={-.1,0,.1}; //this array sets the values that the electric field will take on. ex: -1*X(operator). 
	double e_field_max = 1.0; //upper & lower domain of the electric field. 
	double e_field_min = -1.0;
 
	//double p[3] = {0,0,0}; //initializing the polarization vector. 
	//int m = 3;  //this tells the giant for-loop how many times to go through and corresponds to the number of electric field values. 

	int row = 0; //setting the row value for the chebyshev differentiation matrix. 
	int col = 0; //setting the column value for the chebyshev differentiation matrix. 

	//************** INT M(AND THUS THE NUMBER OF VALUES FOR THE ELECTRIC FIELD) AND NUM_CHEBPOINTS MUST BE THE SAME***********
	double sum; 
	double green;
	double BETA = 0;

	//defining the constant alpha(for the laplacian) and num_chebpoints
	int num_chebpoints = 5; 
	double alpha = (pow(HD -> num_points + 1.0, 2.0)) / (pow(HD -> right_endpoint - HD -> left_endpoint,2.0)); 
	double omega = 1.0;
	int flag = 1.0; 
	int tag = 0;



//******************************************************************************************************

	double *main_diag; //main diagonal for the lapack solver
	main_diag = (double*)malloc(num_points*sizeof(double));

	double *sub_diag; //sub diag for the lapack solver(since this matrix is symetric the sub and super diagonals are the same. 
	sub_diag = (double*)malloc(num_points*sizeof(double));

	double *Z; // identity matrix for the lapacke solver 
	Z = (double*)malloc(num_points*num_points*sizeof(double));

	double *psi; //a place to store the first eigenvector.  
	psi = (double*)malloc(num_points*sizeof(double));

	double *cheb_points;//place to store the chebpoints 
	cheb_points = (double*)malloc(num_chebpoints*sizeof(double)); 

	double *elec_field; //place to store the electric field points which correspond to the chebpoints. 
	elec_field = (double*)malloc(num_chebpoints*sizeof(double)); 

	double *diff; //place to store the value of 
	diff = (double*)malloc(sizeof(double));

	double *polarization; //the polarization vector. 
	polarization = (double*)malloc(num_chebpoints*sizeof(double)); 

	//this assigns the electric field array chebpoints on the user defined e_min & e_max so that was it can vary the electric field using the chebyshev points. 
	cheb_array(num_chebpoints, &elec_field[0], e_field_max, e_field_min); 

	//for(ii<0;ii<num_chebpoints;ii++){
	//	printf("these are the chebpoints: %E \n", elec_field[ii]);
	//}



	//creating giant for loop in order to vary the electric field each time. **************
	for(kk=0;kk<num_chebpoints;kk++){

	//*****creating the hamiltonian*****
 
	lap_explicit(HD -> num_points, &main_diag[0], &sub_diag[0], alpha, 0);

	x_squared(HD -> num_points, HD-> left_endpoint, HD-> right_endpoint, &main_diag[0], omega, 1);

	E_field( HD -> num_points, HD -> left_endpoint, HD-> right_endpoint, &main_diag[0], elec_field[kk] , 1);

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

	LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I', num_points, main_diag, sub_diag, Z, num_points); 


	
	for(ii=0;ii<1;ii++){
	printf("first eigenvalue. %E \n", main_diag[ii]);
	}
	/*
	pick off the first eigenvector e.g. the ground state of the Hamiltonian. 
	for(ii=0;ii<num_points;ii++){
		//psi[ii] = Z[ii]; 




	printf("\nThis is Psi:\n");
	for(ii=0;ii<num_points;ii++){
		printf("Psi[%d] = %+1.15e\n",ii,psi[ii]);
	}
	printf("\n");/*/

	
	normalize_polarize(&Z[0] ,HD -> left_endpoint ,HD-> right_endpoint ,num_points ,&integral);

	polarization[kk] = integral; 

	//printf("The polarization: %1.1e \n", integral);	


}//end of the giant for loop creating the polarization vector.  **************


	for(kk=0;kk<num_chebpoints;kk++){
	printf("The polarization: %E \n", polarization[kk]);
	}


	//creating the array full of cheb_points. 
	//cheb_array(num_chebpoints, &cheb_points[0], e_field_max, e_field_min); 

	/*
	printf("These are the chebpoints: \n");
	for(ii=0;ii<num_chebpoints;ii++){
		printf("cheb_points[ii] = %+1.15e \n", cheb_points[ii]);
	}
	printf("\n");
	/*/
	//creating the first derivative chebyshev differentiation matrix. 
	green = chebdiff(row, col, num_chebpoints, &elec_field[0]); 

 	printf("this is D. %+1.15e \n",green);

	sum = chebdiff2(row, col, num_chebpoints, &elec_field[0]);


 	printf("this is D2 %+1.15e \n",sum);

	BETA = beta(num_chebpoints, &elec_field[0], &polarization[0]); //&elec_field[0] was &cheb_points[0]

	printf("this is beta %1.15e \n", BETA); 
	
	/*
	BETA = 0;

	BETA = beta(row, col, num_chebpoints, &elec_field[0], &polarization[0]);

	printf("This is beta %E \n", BETA); 
	/*/
	free(main_diag);
	free(sub_diag);
	free(cheb_points); 
	free(diff); 

	return 0; 
} 



//*************************************************************************************************
/*
void lap_explicit(int num_points, double *main_diag, double *sub_diag, double alpha, int flag){

	int ii; 

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


void x_squared(int num_points, double left_endpoint, double right_endpoint, double *main_diag, double alpha, int flag){
	
	double interval_length = ( right_endpoint -  left_endpoint) / ( num_points + 1.0);
	int ii; 	

	if(flag == 0){
	for(ii=0;ii<num_points;ii++){
		main_diag[ii] = (pow(( left_endpoint + (ii+1.0)*interval_length),2.0)) * alpha;		
	}

	}

	else{
		for(ii=0;ii<num_points;ii++){
		main_diag[ii] += (pow(( left_endpoint + (ii+1.0)*interval_length),2.0)) * alpha;		
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

//this function should give the chebpoints for an arbitrary interval [e_field min, e_field_max] 
void cheb_array(int num_chebpoints, double *cheb_points, double e_field_max, double e_field_min){

	int ii;

	for(ii=0;ii<num_chebpoints;ii++){

	cheb_points[ii] = -1.0*(((e_field_max - e_field_min))/(2.0)) * sin(M_PI*((ii/(num_chebpoints-1.0)) - 0.5)) + ((e_field_min + e_field_max)/(2.0));	
	}

}




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


double chebdiff2(int row, int col, int num_chebpoints, double *cheb_points){


	int middle = num_chebpoints / 2.0;
	int ii;
	int jj = col; 
	double sum = 0; 

	//for(jj=0;jj<num_chebpoints;jj++){
		for(ii=0;ii<num_chebpoints;ii++){

		sum += chebdiff(middle, ii, num_chebpoints, cheb_points) * chebdiff(ii, jj, num_chebpoints, cheb_points); 
	
		}
	
	return sum; 

}


double beta(int row, int col, int num_chebpoints, double *cheb_points, double *integral){

	int ii; 
	double beta = 0.0;
	
	for(ii=0;ii<num_chebpoints;ii++){
		beta += chebdiff2(row, ii, num_chebpoints, cheb_points ) * integral[ii];
	}

	return beta;
}
/*/







