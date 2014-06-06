//creating the laplacian operator 

//flag should be zero here so that way it will overwrite whatever was in Y. 
//alpha is an arbitrary constant 
//n is number of points
//X is the input vector
//Y is the output vector after the operator 
//flag is whether or not to add the next step to Y or write over Y. 

void lapmult(double alpha ,int num_points, double*X,double*Y,int flag){

	int ii;

	if(flag == 0){
	
	//computing Y using finite difference method where matrix has -1, 2, -1 on tri diagonal. 
	for(ii=1;ii< num_points-1;ii++){
		Y[ii] = (-X[ii-1]+2*X[ii]-X[ii+1])*alpha;
	}
	//top row of lapl
	Y[0] = (2*X[0]-X[1])*alpha;
	
	//bottom row of laplacian 
	Y[num_points-1] = (2*X[num_points-1]-X[num_points-2])*alpha;

	
	//multiplying all terms by 1/h^2 which is alpha. alpha is specified in HamOp so that way it can be an arbitrary scalar. 
	for(ii=0;ii< num_points;ii++){
		Y[ii] = Y[ii]*alpha;
	}

	}

	else{//if else do this

	for(ii=1;ii< num_points -1;ii++){
		Y[ii] += (-X[ii-1]+2*X[ii]-X[ii+1])*alpha;
	}
	
	Y[0] += (2*X[0]-X[1])*alpha;

	Y[num_points-1] += (2*X[ num_points -1]-X[num_points-2])*alpha;

	/* the old way of multiplying by a constant 
	for(ii=0;ii< n ;ii++){
		Y[ii] = Y[ii]*alpha;

	}/*/

	}

	
}

