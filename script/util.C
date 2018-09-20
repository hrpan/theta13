#include "Math/Interpolator.h"

double pull(double x, double mean, double stddev){
	double tmp = ( x - mean ) / stddev;
	return tmp * tmp;
}

double relu(double x){
	return x > 0 ? x : 0;
}

double bound_par(double x, double low, double up){
	if(x < low)
		return low;
	else if(x > up)
		return up;
	else 
		return x;
}


