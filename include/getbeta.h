//header file for creating potential and laplacian and then adding them together to create the Hamiltonian. 
#ifndef __getbeta_h__
#define __getbeta_h__


void lapmult(double alpha, int n, double*X, double*Y, int flag);	

void potential(double beta, int n, double a, double b, double*X, double*Y, int flag); 


#endif 
