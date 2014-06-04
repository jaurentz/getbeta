//this function will operate on an arbitrary vector and multiply by x. e.g. this is the 'X' operator. 

void X_operator(double beta, int n, double a, double b, double* psi_hat, double* function_1, int flag){

	int jj;
	double h = (a - b)/(n+1.0);

	if(flag == 0){	
	//jj+1 because it picks the first interval + 1
	for(jj=0;jj<n;jj++){
		psi_hat[jj] = ( b + (jj+1.0)*h) * (function_1[jj])*beta; 
	}
	
	}
	else{  //if else do this
	for(jj=0;jj<n;jj++){
		psi_hat[jj] += ( b +(jj+1.0)*h) * (function_1[jj])*beta;
	}

	}
 

}

