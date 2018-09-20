#include "constants.h"

bool debug = false;

const double reactor_power[6] = {
	2.514, 2.447, 2.566, 2.519, 2.519, 2.550
};

const double osc_par_inits[4][3] = {
	{0.1, 0, 0.2},
	{0.307, 0, 0.6},
	{7.53e-5, 0, 0},
	{2.43e-3, 0, 0},
};

const double minuit_print = 2;
const double minuit_strategy = 2;

const double simplex_args[2] = {1e5, 0.01};
const double seek_args[2] = {1e4, 1.0};
const double minimize_args[2] = {1e6, 1e-12};
const double minos_args[1] = {1e5};

const bool normal_order = true;
const size_t energy_bins = 10;

const double energy_min = 1.8;
const double energy_max = 15.8;
const double energy_step = (energy_max - energy_min) / energy_bins;

const double det_eff[3] = {0.806, 0.0193, 0.0013};

double ibd_rate[8][2];

vector<double> cross_section;

const double cs_min = 1.8;
const double cs_max = 15.8;
const double cs_step = 0.01;


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

	return 1 - _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 - sin2_2theta13 * _sin2_delta_ee;
}

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

void chi2(int &npar, double *g, double &result, double *par, int flag){
	int par_idx = 0;
	//osc. pars
	double _sin2_2theta13 = par[par_idx++];
	double _sin2_theta12 = par[par_idx++];
	double _sin_theta12 = sqrt(_sin2_theta12);
	double _delta_m2_21 = par[par_idx++];
	double _delta_m2_32 = par[par_idx++];
	double _delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;

	//reactor pars
	double slope0 = par[par_idx++]; 
	double slope1 = par[par_idx++]; 
	double _reactor_shape[energy_bins];
	double _tmp = 0;
	for(size_t bin = 0; bin < energy_bins; ++bin){
		double _e = ( 0.5 + bin ) * energy_step;
		double _p = exp( - slope0 * _e - slope1 * _e * _e );
		_tmp += _p;
		_reactor_shape[bin] = _p;
	}

	double _cross_section[energy_bins];
	for(size_t bin = 0; bin < energy_bins; ++bin){
		 _reactor_shape[bin] /= _tmp;
		_cross_section[bin] = cross_section[int(((0.5 + bin) * energy_step + energy_min - cs_min) / cs_step)];
		if(debug)
			cout << _reactor_shape[bin] << " " << _cross_section[bin] << endl;
	}
	

	//detector pars
	double _det_eff_uncorr[8];
	for(int i=0;i<8;++i)
		_det_eff_uncorr[i] = par[par_idx++];
	double _det_eff_corr = par[par_idx++];
	double _norm = par[par_idx++];
	
	double _chi2 = 0;
	double _penalty = 0;
	for(size_t det = 0; det < 8; ++det){
		double expected = 0;
		for(size_t core = 0; core < 6; ++core){
			double _l = dist_map[det][core];
			double _density = 1 / ( _l * _l );
			for(size_t bin = 0; bin < energy_bins; ++bin){
				double _e = ( 0.5 + bin ) * energy_step + energy_min;
				double _prob = survival_prob(_sin_theta12, _sin2_2theta13, _delta_m2_21, _delta_m2_ee, _l, _e);
				double _cs = _cross_section[bin];
				if(debug)
					cout << bin << " " << _e << " " << int((_e - cs_min)/cs_step) << " " << _prob << endl;
				expected += _prob * _density * reactor_power[core] * _reactor_shape[bin] * _norm * _cross_section[bin];
			}
		}
		expected *= det_eff[0] * (1 + _det_eff_corr + _det_eff_uncorr[det]);
		if(debug)
			cout << expected << " " << ibd_rate[det][0] << " " << pull(expected, ibd_rate[det][0], ibd_rate[det][1]) << endl;
		_chi2 += pull(expected, ibd_rate[det][0], ibd_rate[det][1]);

		_penalty += pull(_det_eff_uncorr[det], 0, det_eff[2]);
	}		
	_penalty += pull(_det_eff_corr, 0, det_eff[1]);
	_penalty += pull(_sin2_theta12, sin2_theta12[0], sin2_theta12[1]);
	_penalty += pull(_delta_m2_21, delta_m2_21[0], delta_m2_21[1]);
	_penalty += pull(_delta_m2_32, delta_m2_32[0], delta_m2_32[1]);
	result = _chi2 + _penalty;
}

void parseInput(){
	double _nus[8];
	for(int i=0;i<8;++i) cin >> _nus[i];
	double _lt[8];
	for(int i=0;i<8;++i) cin >> _lt[i];
	double _eps_mu[8];
	for(int i=0;i<8;++i) cin >> _eps_mu[i];
	double _eps_mult[8];
	for(int i=0;i<8;++i) cin >> _eps_mult[i];
	double _bkgs[6][8][2];
	for(int b=0;b<6;++b) for(int i=0;i<8;++i) cin >> _bkgs[b][i][0] >> _bkgs[b][i][1];

	for(int i=0;i<8;++i) cin >> ibd_rate[i][0] >> ibd_rate[i][1];
}

void parseCrossSection(){
	ifstream ifile("./data/cross_section");
	double _e, _cs;
	while( ifile >> _e >> _cs ){
		cross_section.push_back(_cs * 1e44);
	}
}

void fit(){
	parseInput();		
	parseCrossSection();

	TFitter minuit(16);
	minuit.ExecuteCommand("SET PRINTOUT", &minuit_print, 1);
	minuit.ExecuteCommand("SET STRATEGY", &minuit_strategy, 1);
	minuit.SetFCN(chi2);

	int idx = 0;
	for(int p=0;p<4;++p){
		char buf[255];
		sprintf(buf, "OSC_PAR%d", p);
		minuit.SetParameter(idx++, buf, osc_par_inits[p][0], 0.01 * osc_par_inits[p][0], osc_par_inits[p][1], osc_par_inits[p][2]);
	}
	for(int p=0;p<2;++p){
		char buf[255];
		sprintf(buf, "REACTOR_SHAPE%d", p);
		minuit.SetParameter(idx++, buf, 1, 1e-3, 0, 2);
	}
	minuit.SetParameter(5, "TEST", 0, 0, 0, 0);
	minuit.FixParameter(5);
	for(int p=0;p<8;++p){
		char buf[255];
		sprintf(buf, "DET_EFF_UNCORR%d", p);
		minuit.SetParameter(idx++, buf, 0, 1e-4, 0, 0);
	}
	minuit.SetParameter(idx++, "DET_EFF_CORR", 0, 1e-4, 0, 0);
	minuit.SetParameter(idx++, "NORM", 1e4, 1e2, 0, 0);
	if(!debug){
		minuit.ExecuteCommand("SIMPLEX", simplex_args, 2);
		//minuit.ExecuteCommand("SEEK", seek_args, 2);
		minuit.ExecuteCommand("MINIMIZE", minimize_args, 2);
	}else{
		double __tmp = 1;
		minuit.ExecuteCommand("SEEK", &__tmp, 1);
	}
//	minuit.ExecuteCommand("MINOS", minos_args, 1);
}
