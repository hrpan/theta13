#include "Math/Interpolator.h"

void parseInput(double ibd_rate[8][2]){
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

void parseCrossSection(vector<double> &cross_section, double &e_min, double &e_max, double &e_step){
	ifstream ifile("./data/cross_section");
	double _e, _cs;
	while( ifile >> _e >> _cs ){
		e_min = min(_e, e_min);
		e_max = max(_e, e_max);
		cross_section.push_back(_cs * 1e44);
	}
	e_step = ( e_max - e_min ) / cross_section.size();
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
