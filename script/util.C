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

double d_ji(double d_m2_ji, double l, double e_nu){
	return 1.267 * d_m2_ji * l / e_nu;
}

double survival_prob(double sin_theta12, double sin2_2theta13, double d_21, double d_ee, double l, double e_nu){
	double _theta13 = asin(sqrt(sin2_2theta13)) / 2;
	double _cos_theta13 = cos(_theta13);
	double _cos4_theta13 = _cos_theta13 * _cos_theta13 * _cos_theta13 * _cos_theta13;
	double _theta12 = asin(sin_theta12); 
	double _sin2_2theta12 = sin(2 * _theta12) * sin(2 * _theta12);
	double _sin2_delta_21 = sin(d_ji(d_21, l, e_nu)) * sin(d_ji(d_21, l, e_nu));
	double _sin2_delta_ee = sin(d_ji(d_ee, l, e_nu)) * sin(d_ji(d_ee, l, e_nu));
	//cout << _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 << " " << sin2_2theta13 * _sin2_delta_ee << " " << _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 / (sin2_2theta13 * _sin2_delta_ee) << endl;
	return 1 - _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 - sin2_2theta13 * _sin2_delta_ee;
}

TGraph *minuit_profile(TFitter &minuit, int par_no, int npoints, double ndevs){
	cout << "STARTING TO PROFILE PARAMETER" << par_no << endl;
	gMinuit->SetPrintLevel(-1);

	double _p0 = minuit.GetParameter(par_no);
	double _p0_dev = minuit.GetParError(par_no);	
	
	double amin, edm, errdef;
	int nvpar, nparx;
	minuit.GetStats(amin, edm, errdef, nvpar, nparx);

	double chi2_min = amin;

	minuit.FixParameter(par_no);
	vector<double> x, y;
	double scan_step = scan_dev / scan_pts;
	for(int pt = - (npoints - 1); pt < npoints; ++pt){
		cout << "\nSCANNING PAR" << par_no << " AT " << pt * scan_step << "DEV" << endl;
		double _tmp[2] = {par_no + 1, _p0 + _p0_dev * scan_step * pt};
		minuit.ExecuteCommand("SET PARAMETER", _tmp, 2);
		minuit.ExecuteCommand("MINIMIZE", 0, 0);
		minuit.GetStats(amin, edm, errdef, nvpar, nparx);
		x.push_back(_tmp[1]);
		y.push_back(amin - chi2_min);
	}
	minuit.ReleaseParameter(par_no);	

	TGraph *g = new TGraph(x.size(), &x[0], &y[0]);
	return g;
}
