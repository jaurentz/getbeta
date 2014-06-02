//creating the laplacian operator 

void lapmult(double alpha ,int n, double*X,double*Y,int flag){

	int ii;

	if(flag == 0){
	
	//computing Y
	for(ii=1;ii< n-1;ii++){
		Y[ii] = (-X[ii-1]+2*X[ii]-X[ii+1])*alpha;
	}
	//top row of lapl
	Y[0] = (2*X[0]-X[1])*alpha;
	
	//bottom row of laplacian 
	Y[n-1] = (2*X[n-1]-X[n-2])*alpha;

	
	//multiplying all terms by 1/h^2
	for(ii=0;ii< n;ii++){
		Y[ii] = Y[ii]*alpha;
	}

	}

	else{//if else do this

	for(ii=1;ii< n-1;ii++){
		Y[ii] += (-X[ii-1]+2*X[ii]-X[ii+1])*alpha;
	}
	
	Y[0] += (2*X[0]-X[1])*alpha;

	Y[n-1] += (2*X[ n-1]-X[n-2])*alpha;

	/* the old way of multiplying by a constant 
	for(ii=0;ii< n ;ii++){
		Y[ii] = Y[ii]*alpha;

	}/*/

	}

	
}

