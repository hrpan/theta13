#include "util.C"
#include "constants.h"

double chi2_stats(
	double _sin2_2theta13,
	double _sin_theta12,
	double _delta_m2_21,
	double _delta_m2_ee,
	double *_reactor_shape,
	double _weighted_reactor_power[n_det][n_core],
	double *_cross_section,
	size_t energy_bins,
	double *_det_eff_uncorr,
	double _det_eff_corr,
	double *_bkg_uncorr,
	double _norm,
	double obs_rate[n_det][2],
	bool debug,
	bool print
	){

	
	double _chi2 = 0;
	for(size_t det = 0; det < n_det; ++det){
		double expected = 0;
		for(size_t core = 0; core < n_core; ++core){
			double _l = dist_map[det][core];
			double _density = 1 / ( _l * _l );
			double _sum = 0;
			for(size_t bin = 0; bin < energy_bins; ++bin){
				double _e = ( 0.5 + bin ) * energy_step + energy_min;
				double _prob = survival_prob(_sin_theta12, _sin2_2theta13, _delta_m2_21, _delta_m2_ee, _l, _e);
				double _n_bin = _prob * _reactor_shape[bin] * _cross_section[bin];
				if(debug)
					cout << "bin:" << bin << " e:" << _e << " prob:" << _prob << " r_shape:" << _reactor_shape[bin] << " c_sect:" << _cross_section[bin] << " n_bin:" << _n_bin << endl;
				_sum += _n_bin;
				expected += _prob * _density * _weighted_reactor_power[det][core] * _reactor_shape[bin] * _norm * _cross_section[bin] * energy_step;
			}
			_sum *= (energy_step * _density * _norm);

			if(debug)
				cout << "sum:" << _sum << endl;

			expected += _sum; 
		}

		expected *= (det_eff[0] + _det_eff_corr + _det_eff_uncorr[det]);

		if(debug)
			cout << "exp:"<< expected << " obs:" << obs_rate[det][0] << " pull:" << pull(expected + _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1]) << endl;

		if(print){
			printf("det: %d exp: %10.5f bkg: %10.5f obs: %10.5f/%10.5f chi2: %10.5f\n", det, expected, _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1], pull(expected + _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1]));
			fflush(stdout);
		}
		_chi2 += pull(expected + _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1]);

	}	
	return _chi2;	
}

double chi2_syst(
	double _sin2_theta12,
	double _delta_m2_21,
	double _delta_m2_32,
	double *_det_eff_uncorr,
	double _det_eff_corr,
	double *_bkg_uncorr,
	bool print
	){
	
	double _chi2 = 0;
	for(size_t det = 0; det < 8; ++det){
		_chi2 += pull(_bkg_uncorr[det], bkg_rate[det][0], bkg_rate[det][1]);
		_chi2 += pull(_det_eff_uncorr[det], 0, det_eff[2]);
		if(print)
			printf("det: %d bkg_uncorr: %10.5f det_eff_uncorr: %10.5f\n", det, pull(_bkg_uncorr[det], bkg_rate[det][0], bkg_rate[det][1]), pull(_det_eff_uncorr[det], 0, det_eff[2]));
	}
	_chi2 += pull(_det_eff_corr, 0, det_eff[1]);
	_chi2 += pull(_sin2_theta12, sin2_theta12[0], sin2_theta12[1]);
	_chi2 += pull(_delta_m2_21, delta_m2_21[0], delta_m2_21[1]);
	_chi2 += pull(_delta_m2_32, delta_m2_32[0], delta_m2_32[1]);
	if(print){
		printf("det_eff_corr: %10.5f\n", pull(_det_eff_corr, 0, det_eff[1]));
		printf("sin2_theta12: %10.5f\n", pull(_sin2_theta12, sin2_theta12[0], sin2_theta12[1]));
		printf("delta_m2_21: %10.5f\n", pull(_delta_m2_21, delta_m2_21[0], delta_m2_21[1]));
		printf("delta_m2_32: %10.5f\n", pull(_delta_m2_32, delta_m2_32[0], delta_m2_32[1]));
		fflush(stdout);
	}
	return _chi2;	
}


