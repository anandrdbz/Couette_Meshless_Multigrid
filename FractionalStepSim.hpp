#ifndef FRAC_STEP_SIM_H
#define FRAC_STEP_SIM_H
#include "fractionalStepGrid.hpp"
#include "testing_functions.hpp"
#include "FracStepMultigrid.hpp"
class FractionalStepParams {
public:
	std::string flowType;
	double rho, dt, mu, ppe_conv_res;
	vector<string> filenames;
	vector<GridProperties> props;
	string extension;
	string directory;

};
FractionalStepParams gen_fracstep_param(int numGrids, int poly_deg, double dt, double mu, double rho, double ppe_conv);
FractionalStepGrid* genFractionalStepGrid(const char* filename, GridProperties props, double dt, double mu, double rho, double ppe_conv, std::string coarse, bool gmres, bool regularize);
void check_derivs(FractionalStepGrid* grid);
void run_fracstep_param(FractionalStepParams params, double endtime, std::string Solver);
void run_frac_step_test();

#endif