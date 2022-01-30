#ifndef FRAC_STEP_GRID_H
#define FRAC_STEP_GRID_H
#include "grid.h"
class FractionalStepGrid : public Grid {
public:

	double dt;
	double ppe_conv_res;
	double rho;
	double mu;
	double lambda;
	std::string flowType;
	Eigen::VectorXd *u, *v, *u_old, *v_old, *v_hat, *u_hat;
	Eigen::SparseMatrix<double, Eigen::RowMajor> *derivXMat_, *derivYMat_, *uvLaplaceMat_,*derivXMat2_,*derivYMat2_,*uvLaplaceMat2_;

	FractionalStepGrid(vector<std::tuple<double, double, double>> points, vector<Boundary> boundaries,
		GridProperties properties, Eigen::VectorXd source);
	~FractionalStepGrid();
	void set_uv_bound();
	void build_derivX_mat();
	void build_derivY_mat();
	void build_uv_laplace_mat();
	void calc_u_hat();
	void calc_v_hat();
	void set_ppe_source();
	void correct_u();
	void correct_v();
	double fs_residual();
	double delu();
	void prescribe_soln();
};

#endif