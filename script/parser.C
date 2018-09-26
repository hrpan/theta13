#include "Math/Interpolator.h"

void parseInput(double ibd_rate[8][2], double bkg_rate[8][2], double obs_rate[8][2]){
	double _nus[8];
	for(int i=0;i<8;++i) cin >> _nus[i];
	double _lt[8];
	for(int i=0;i<8;++i) cin >> _lt[i];
	double _eps_mu[8];
	for(int i=0;i<8;++i) cin >> _eps_mu[i];
	double _eps_mult[8];
	for(int i=0;i<8;++i) cin >> _eps_mult[i];
	double _bkgs[6][8][2];
	for(int b=0;b<6;++b) 
		for(int i=0;i<8;++i){ 
			cin >> _bkgs[b][i][0] >> _bkgs[b][i][1];
			bkg_rate[i][0] += _bkgs[b][i][0];
			bkg_rate[i][1] += ( _bkgs[b][i][1] * _bkgs[b][i][1] );
		}
	for(int i=0;i<8;++i){
		bkg_rate[i][1] = sqrt(bkg_rate[i][1]);
		cout << bkg_rate[i][0] << " " << bkg_rate[i][1] << " ";
		double _tmp = _lt[i] * _eps_mu[i] * _eps_mult[i];
		obs_rate[i][0] = _nus[i]/ _tmp;
		obs_rate[i][1] = sqrt(_nus[i]) / _tmp;
		cout << obs_rate[i][0] << " " << obs_rate[i][1] << " ";
	}
	for(int i=0;i<8;++i) cin >> ibd_rate[i][0] >> ibd_rate[i][1];
}

void parseCrossSection(ROOT::Math::Interpolator &*f, double &e_min, double &e_max){
	ifstream ifile("./data/cross_section");
	double _e, _cs;
	vector<double> x, y;
	while( ifile >> _e >> _cs ){
		x.push_back(_e);
		y.push_back(_cs * 1e44);
		e_min = min(_e, e_min);
		e_max = max(_e, e_max);
	}
	f = new ROOT::Math::Interpolator(x, y, 0);
}

void parseReactorSpectrum(const char *data, ROOT::Math::Interpolator &*f){
	ifstream ifile(data);
	string str;
	vector<double> x, y;
	while(getline(ifile, str)){
		if(str[0] == '#') 
			continue;

		istringstream tmp(str);
		
		double _x, _y;
		tmp >> _x >> _y;
	//	cout << _x << " " << _y << endl;			
		x.push_back(_x);
		y.push_back(_y);
	}
	f = new ROOT::Math::Interpolator(x, y, 0);
}
