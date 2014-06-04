//this product will take two different vectors and compute the inner product over the interval. This is the same as the integral of psi^2 but this function can take two different vectors or the same vector. 

//This function relies on Trapz function 

void inner_product(double * function_1, double * function_2, double left_endpoint, double right_endpoint, int interior_points, double * integral){
	
	int ii; 

	double * psi_hat;
	psi_hat = (double*)malloc(interior_points*sizeof(double));

	for(ii=0;ii<interior_points;ii++){
		psi_hat[ii] = function_1[ii] * function_2[ii];
	}

	trapz(psi_hat, left_endpoint, right_endpoint, interior_points, integral);

	free(psi_hat);

}

