//header file for creating potential and laplacian and then adding them together to create the Hamiltonian. 
#ifndef __getbeta_h__
#define __getbeta_h__


void lapmult(double alpha, int n, double*X, double*Y, int flag);	

void potential(double beta, int n, double a, double b, double*X, double*Y, int flag); 

//these are the integration functions

void trapz(double * function_values, double left_endpoint, double right_endpoint, int interior_points, double *integral);   

void inner_product(double * function_1, double * function_2, double left_endpoint, double right_endpoint, int interior_points, double *integral);

void polarization(double * function_1,double left_endpoint, double right_endpoint, int interior_points, double * integral); 

void X_operator(double beta, int n, double a, double b, double *psi_hat, double *function_1, int flag);

void norm(double * function_1, double left_endpoint, double right_endpoint, int interior_points, double * integral);



#endif 
