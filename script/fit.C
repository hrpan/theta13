#include "constants.h"
#include "parser.C"
#include "util.C"

bool debug = false;

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
const double minimize_args[2] = {1e6, 1e-8};
const double minos_args[1] = {1e5};

const size_t energy_bins = 50;

const double energy_min = 2;
const double energy_max = 8;
const double energy_step = (energy_max - energy_min) / energy_bins;

const double det_eff[3] = {0.806, 0.0193, 0.0013};

ROOT::Math::Interpolator *r_spec[reactor_nuclears];

double ibd_rate[8][2];

vector<double> cross_section;
double cs_min, cs_max, cs_step;


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
	double _r_frac[reactor_nuclears];
	double _r_frac_sum = 1;
	_r_frac[reactor_nuclears-1] = 1;
	for(size_t n = 0; n < reactor_nuclears-1; ++n){
		double _f = par[par_idx++];
		_r_frac[n] = _f;
		_r_frac_sum += _f;
	}
	for(size_t n = 0; n < reactor_nuclears; ++n){
		_r_frac[n] /= _r_frac_sum;
		if(debug)
			cout << "_r_frac" << n << " " << _r_frac[n] << endl;
	}
	double _reactor_shape[energy_bins];
	double _tmp = 0;
	for(size_t bin = 0; bin < energy_bins; ++bin){
		double _e = ( 0.5 + bin ) * energy_step + energy_min;
		double _sum = 0;
		for(size_t n = 0; n < reactor_nuclears; ++n)
			_sum += _r_frac[n] * r_spec[n]->Eval(_e);
		_tmp += _sum;
		_reactor_shape[bin] = _sum;
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
				if(debug)
					cout << "bin" << bin << " " << _e << " " << _prob << " " << _density << " " << _reactor_shape[bin] << endl;
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

void fit(){
	parseInput(ibd_rate);		
	parseCrossSection(cross_section, cs_min, cs_max, cs_step);

	for(size_t r = 0; r < reactor_nuclears; ++r)
		parseReactorSpectrum(reactor_data[r], r_spec[r]);

	TFitter minuit(14 + reactor_nuclears - 1);
	minuit.ExecuteCommand("SET PRINTOUT", &minuit_print, 1);
	minuit.ExecuteCommand("SET STRATEGY", &minuit_strategy, 1);
	minuit.SetFCN(chi2);

	int idx = 0;
	char buf[255];
	for(int p=0;p<4;++p){
		sprintf(buf, "OSC_PAR%d", p);
		minuit.SetParameter(idx++, buf, osc_par_inits[p][0], 0.01 * osc_par_inits[p][0], osc_par_inits[p][1], osc_par_inits[p][2]);
	}
	for(int p=0;p<reactor_nuclears-1;++p){
		sprintf(buf, "REACTOR_SHAPE%d", p);
		minuit.SetParameter(idx++, buf, 1, 0.1, 0, 100);
	}
	for(int p=0;p<8;++p){
		sprintf(buf, "DET_EFF_UNCORR%d", p);
		minuit.SetParameter(idx++, buf, 0, 1e-4, 0, 0);
	}
	minuit.SetParameter(idx++, "DET_EFF_CORR", 0, 1e-4, 0, 0);
	minuit.SetParameter(idx++, "NORM", 1e5, 1e2, 0, 0);
	if(!debug){
		minuit.ExecuteCommand("SIMPLEX", simplex_args, 2);
		//minuit.ExecuteCommand("SEEK", seek_args, 2);
		minuit.ExecuteCommand("MINIMIZE", minimize_args, 2);
	}else{
		double __tmp = 5;
		minuit.ExecuteCommand("SEEK", &__tmp, 1);
	}
//	minuit.ExecuteCommand("MINOS", minos_args, 1);
}
