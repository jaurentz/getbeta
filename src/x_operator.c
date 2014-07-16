//this function will operate on an arbitrary vector and multiply by x. e.g. this is the 'X' operator. 


void X_operator(double beta, int num_points, double left_endpoint, double right_endpoint, double* psi_hat, double* function_1, int flag){


	int jj;
	double interval_length = (right_endpoint - left_endpoint)/(num_points + 1.0);

	if(flag == 0){	
	//jj+1 because it picks the first interval + 1
	for(jj=0;jj<num_points;jj++){
		psi_hat[jj] = ( left_endpoint + (jj+1.0)* interval_length ) * (function_1[jj])*beta; 
	}
	
	}
	else{  //if else do this
	for(jj=0;jj<num_points;jj++){
		psi_hat[jj] += ( left_endpoint +(jj+1.0)* interval_length) * (function_1[jj])*beta;
	}

	}
 

}

