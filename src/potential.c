#include <math.h>

//this function creates(computes) the potential operator for an arbitary X vector(its address) and outputs a Y(its address) vector.  

//flag should be 1 here because this should be added to the value of Y now because lapmult will have already acted on Y. So, you do not want to overwrite Y. 

//computing h, the size of the interval, locally in this function. 

//beta is computed in HamOp function and is an arbitrary scalar. 



void potential(double beta, int num_points, double right_endpoint, double left_endpoint, double*X, double*Y,int flag){

	int jj;
	double interval_length = (right_endpoint - left_endpoint)/(num_points + 1.0);

	if(flag == 0){	
	//jj+1 because it picks the first interval + 1
	for(jj=0;jj<num_points;jj++){
		Y[jj] = (pow(( left_endpoint + (jj+1.0)*interval_length),2) * X[jj])*beta; 
	}
	
	} //note that here Y[jj] is plus/equals which means that it will add to the previous value of Y[] instead of overwriting with its value.
	else{  //if else do this
	for(jj=0;jj<num_points;jj++){
		Y[jj] += (pow(( left_endpoint +(jj+1.0)*interval_length),2) * X[jj])*beta;
	}

	}
 

}

