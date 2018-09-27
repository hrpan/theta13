#include "constants.h"
#include "parser.C"
#include "util.C"
#include "chi2.C"

bool debug = false;
bool draw_contour = false;


const double osc_par_inits[4][3] = {
	{8.60e-2, 0, 0.2},
	{3.07e-1, 0, 1.0},
	{7.53e-5, 0, 0},
	{2.43e-3, 0, 0}
};

const size_t npars = 24;

const int scan_pts = 30;
const double scan_dev = 5;

const double minuit_print = 1;
const double minuit_strategy = 1;

const double simplex_args[2] = {1e5, 0.01};
const double seek_args[2] = {1e4, 1.0};
const double minimize_args[2] = {1e6, 1e-3};
const double minos_args[1] = {1e5};

const size_t energy_bins = 10;

const double energy_min = 2;
const double energy_max = 8;
const double energy_step = (energy_max - energy_min) / energy_bins;

ROOT::Math::Interpolator *r_spec[reactor_nuclears];
ROOT::Math::Interpolator *cross_section;

double _cross_section[energy_bins];

double _weighted_reactor_power[8][6];

double ibd_rate[8][2], bkg_rate[8][2], obs_rate[8][2];

double cs_min, cs_max;

void parseParameters(
	double *par,
	double &_sin2_2theta13, 
	double &_sin2_theta12,
	double &_sin_theta12, 
	double &_delta_m2_21,
	double &_delta_m2_32, 
	double &_delta_m2_ee,
	double _reactor_shape[energy_bins],
	double _det_eff_uncorr[8],
	double &_det_eff_corr,
	double _bkg_uncorr[8],
	double &_norm){

	int par_idx = 0;
	//osc. pars
	_sin2_2theta13 = par[par_idx++];
	_sin2_theta12 = par[par_idx++];
	_sin_theta12 = sqrt(_sin2_theta12);
	_delta_m2_21 = par[par_idx++];
	_delta_m2_32 = par[par_idx++];
	_delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;

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

	for(size_t bin = 0; bin < energy_bins; ++bin){
		double _e = ( 0.5 + bin ) * energy_step + energy_min;
		double _sum = 0;
		for(size_t n = 0; n < reactor_nuclears; ++n)
			_sum += _r_frac[n] * r_spec[n]->Eval(_e);
		_reactor_shape[bin] = _sum;
	}


	//detector pars
	for(int i=0;i<8;++i)
		_det_eff_uncorr[i] = par[par_idx++];
	_det_eff_corr = par[par_idx++];


	//bkg pars
	for(int i=0;i<8;++i)
		_bkg_uncorr[i] = par[par_idx++];

	_norm = par[par_idx++];

}

void chi2_wrap(int &npar, double *g, double &result, double *par, int flag){

	result = chi2(par, false);
			
	printf("chi2: %20.10f\r", result);
	fflush(stdout);
}

double chi2(double *par, bool print){
	double _sin2_2theta13;
	double _sin2_theta12;
	double _sin_theta12;
	double _delta_m2_21;
	double _delta_m2_32;
	double _delta_m2_ee;
	double _reactor_shape[energy_bins];
	double _det_eff_uncorr[8];
	double _det_eff_corr;
	double _bkg_uncorr[8];
	double _norm;
	
	parseParameters(
	 par,
	 _sin2_2theta13, 
	 _sin2_theta12,
	 _sin_theta12, 
	 _delta_m2_21,
	 _delta_m2_32, 
	 _delta_m2_ee,
	 _reactor_shape,
	 _det_eff_uncorr,
	 _det_eff_corr,
	 _bkg_uncorr,
	 _norm);


	double _chi2_stats = chi2_stats(
		_sin2_2theta13, 
		_sin_theta12, 
		_delta_m2_21,
		_delta_m2_ee, 
		_reactor_shape, 
		_weighted_reactor_power, 
		_cross_section, 
		energy_bins, 
		_det_eff_uncorr, 
		_det_eff_corr,
		_bkg_uncorr, 
		_norm, 
		obs_rate, 
		debug, 
		print);

	double _chi2_syst = chi2_syst(
		_sin2_theta12, 
		_delta_m2_21, 
		_delta_m2_32,
		_det_eff_uncorr, 
		_det_eff_corr, 
		_bkg_uncorr, 
		print);

	if(print){
		printf("chi2_stats: %10.5f chi2_syst: %10.5f chi2_sum: %10.5f\n", _chi2_stats, _chi2_syst, _chi2_stats + _chi2_syst);
		fflush(stdout);
	}
	
	return _chi2_stats + _chi2_syst;
}


void fit(){
	parseInput(ibd_rate, bkg_rate, obs_rate);		
	parseCrossSection(cross_section, cs_min, cs_max);
	for(size_t bin = 0; bin < energy_bins; ++bin){
		double _e = ( 0.5 + bin ) * energy_step + energy_min;
		_cross_section[bin] = cross_section->Eval(_e);
	}

	for(size_t r = 0; r < reactor_nuclears; ++r)
		parseReactorSpectrum(reactor_data[r], r_spec[r]);

	for(size_t det = 0; det < 8; ++det){
		for(size_t core = 0; core < 6; ++core){
			if(det == 3 || det == 7)
				_weighted_reactor_power[det][core] = reactor_power_8ad[core];
			else
				_weighted_reactor_power[det][core] = lt_6ad * reactor_power_6ad[core] + (1- lt_6ad) * reactor_power_8ad[core]; 
		}
	}


	TFitter minuit(npars);
	minuit.ExecuteCommand("SET PRINTOUT", &minuit_print, 1);
	minuit.ExecuteCommand("SET STRATEGY", &minuit_strategy, 1);
	minuit.SetFCN(chi2_wrap);

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
	minuit.FixParameter(4);
	minuit.FixParameter(5);
	for(int p=0;p<8;++p){
		sprintf(buf, "DET_EFF_UNCORR%d", p);
		minuit.SetParameter(idx++, buf, 0, 1e-4, 0, 0);
	}
	minuit.SetParameter(idx++, "DET_EFF_CORR", 0, 1e-4, 0, 0);
	for(int p=0;p<8;++p){
		sprintf(buf, "BKG_UNCORR%d", p);
		minuit.SetParameter(idx++, buf, bkg_rate[p][0], bkg_rate[p][1], 0, 0);
	}
	minuit.SetParameter(idx++, "NORM", 2.22e5, 1e3, 0, 0);
	if(!debug){
		//minuit.ExecuteCommand("SIMPLEX", simplex_args, 2);
		//minuit.ExecuteCommand("SEEK", seek_args, 2);
		minuit.ExecuteCommand("MINIMIZE", minimize_args, 2);
	}else{
		double __tmp = 1;
		minuit.ExecuteCommand("SEEK", &__tmp, 1);
	}
	if(draw_contour){
		TGraph *g = (TGraph*)gMinuit->Contour(100, 0, 3);
		g->Draw();
		c1->SaveAs("test.png");
	}
	double pars[npars][2];
	for(size_t n=0;n<npars;++n){
		pars[n][0] = minuit.GetParameter(n);
		pars[n][1] = minuit.GetParError(n);
	}

	TGraph *g = minuit_profile(minuit, 0, scan_pts, scan_dev);
	g->SetMinimum(0);
	g->Draw();

	for(int dev=1; dev < 5; ++dev){
		int line_style = 9;
		double _c = pars[0][0];
		double _d = dev * pars[0][1];
		double _dev2 = dev * dev;
		TLine *_x1 = new TLine(_c - _d, 0, _c - _d, _dev2);
		_x1->SetLineStyle(line_style);
		TLine *_x2 = new TLine(_c + _d, 0, _c + _d, _dev2);
		_x2->SetLineStyle(line_style);
		TLine *_y1 = new TLine(_c - _d, _dev2, _c + _d, _dev2);
		_y1->SetLineStyle(line_style);
		_x1->Draw("same");
		_x2->Draw("same");
		_y1->Draw("same");	
	}

	c1->SaveAs("profile.png");
}
