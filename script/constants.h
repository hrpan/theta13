const size_t n_det = 8;
const size_t n_core = 6;

const double dist_map[n_det][n_core] = {
	{362.38, 371.76, 903.47, 817.16, 1353.62, 1265.32},
	{357.94, 368.41, 903.35, 816.90, 1354.23, 1265.89},
	{1332.48, 1358.15, 467.57, 489.58, 557.58, 499.21},
	{1337.43, 1362.88, 472.97, 495.35, 558.71, 501.07},
	{1919.63, 1894.34, 1533.18, 1533.63, 1551.38, 1524.94},
	{1917.52, 1891.98, 1534.92, 1535.03, 1554.77, 1528.05},
	{1925.26, 1899.86, 1538.93, 1539.47, 1556.34, 1530.08},
	{1923.15, 1897.51, 1540.67, 1540.87, 1559.72, 1533.18},
};

const double reactor_power_6ad[n_core] = {
	2.082, 2.874, 2.516, 2.554, 2.825, 1.976
};

const double reactor_power_8ad[n_core] = {
	2.514, 2.447, 2.566, 2.519, 2.519, 2.550
};

const double lt_6ad = 0.1704238;

const bool normal_order = true;

const double sin2_theta12[2] = {0.307, 0.013};
const double delta_m2_21[2] = {7.53e-5, 0.18e-5}; 
const double delta_m2_32[2] = {2.43e-3, 0.07e-3}; 

const size_t reactor_nuclears = 3;

const char *reactor_data[reactor_nuclears] = {
	"./data/Pu239.dat",
	"./data/Pu241.dat",
	"./data/U235.dat"
};

const double det_eff[3] = {0.806, 0.0193, 0.0013};
