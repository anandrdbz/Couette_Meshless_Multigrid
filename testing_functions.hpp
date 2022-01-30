#ifndef TESTING_FUNCTIONS_H
#define TESTING_FUNCTIONS_H
#include "multigrid.h"
#include "fileReadingFunctions.h"
#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
#define pi 3.141592653589793238462643383279
using std::cout;
using std::endl;
using std::string;

class MultigridParameters {
public:
	vector<string> filenames;
	vector<string> filetypes;
	vector<GridProperties> props;
	string extension;
	string directory;
	string geomtype;
	bool neumann;
	int k1;
	int k2;
	int num_v_cycle;
};
double calc_l1_error(Grid* grid, bool neumannFlag, int k1, int k2);
double calc_l1_error_circle(Grid* grid, bool neumannFlag, int k);
Grid* genGmshGridDirichlet(string geomtype, const char* filename, GridProperties props, std::string filetype, int k1, int k2);
Grid* genGmshGridNeumann(string geomtype, const char* filename, GridProperties props, std::string filetype, int k1, int k2, std::string coarse);
void write_temp_contour (Grid* grid, string directory, std::string extension);
void write_mg_resid(Multigrid& mg, string directory, string extension);
void write_l1error_cond(Multigrid& mg, MultigridParameters params, string directory, string extension, double time);
void run_mg_sim(MultigridParameters params);
MultigridParameters gen_mg_param(string geom, int numGrids, int k, int poly_deg, int vcyc, bool neumann);
void run_tests();
void testGmshSingleGrid();
#endif
