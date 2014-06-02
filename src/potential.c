//creating the potential operator 
void potential(double beta, int n, double a, double b, double*X, double*Y,int flag){

	int jj;
	double h = (a - b)/(n+1.0);

	if(flag == 0){	
	
	for(jj=0;jj<n;jj++){
		Y[jj] = (pow(( b + (jj+1.0)*h),2) * X[jj])*beta; 
	}
	
	}
	else{  //if else do this
	for(jj=0;jj<n;jj++){
		Y[jj] += (pow(( b +(jj+1.0)*h),2) * X[jj])*beta;
	}

	}
 

}

