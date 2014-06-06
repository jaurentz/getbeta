/*header file for creating potential and laplacian and then adding them together to create the Hamiltonian. 
**These functions are for 1-Dimension Schrodinger Equation using finite difference method, and are Homogenous Dirichlet with Boundary conditions that they are zero. This is abbreviated with the term "1DFDHDBC".
/*/
#ifndef __getbeta_h__
#define __getbeta_h__


void lapmult(double alpha, int num_points, double*X, double*Y, int flag);	

void potential(double beta, int num_points, double right_endpoint, double left_endpoint, double*X, double*Y, int flag); 

//these are the integration functions

void trapz(double * function_values, double left_endpoint, double right_endpoint, int num_points, double *integral);   

void inner_product(double * function_1, double * function_2, double left_endpoint, double right_endpoint, int interior_points, double *integral);

void polarization(double * function_1,double left_endpoint, double right_endpoint, int interior_points, double * integral); 

void X_operator(double beta, int num_points, double left_endpoint, double right_endpoint, double *psi_hat, double *function_1, int flag);

void norm(double * function_1, double left_endpoint, double right_endpoint, int interior_points, double * integral);

void normalize_polarize(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral);



#endif 
