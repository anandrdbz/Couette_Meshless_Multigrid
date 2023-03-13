#ifndef GRID_H
#define GRID_H
#include "gridclasses.hpp"
#include <tuple>
#include <string>
#include "fileReadingFunctions.h"
#include "general_computation_functions.h"
#include <stdexcept>
#include "math.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <time.h>

typedef std::tuple<double, double, double> Point;
using std::vector;
using std::cout;
using std::endl;
class Grid {
public:
	//Instance Variables
	Eigen::VectorXd* values_;
	Eigen::VectorXd* sol_;
	Eigen::VectorXd* residuals_;
	Eigen::VectorXd source_;
	Eigen::VectorXd rhs_;
	vector<Point> points_;
	vector<Boundary> boundaries_;
	std::vector<Point> normalVecs_;
	vector<deriv_normal_bc> deriv_normal_coeffs_;
	GridProperties properties_;
	vector<std::pair<int, int>> ptsConn_;
	int laplaceMatSize_;
	int bcount_;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* laplaceMat_;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* restrict_mat;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* prolong_mat;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* interior_mat;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* neumann_boundary_coeffs_;
	Eigen::IncompleteLUT<double,int> solver;
	Eigen::VectorXd diags;
	vector<int> bcFlags_;
	bool neumannFlag_;
	bool regularizeFlag_;
	bool implicitFlag_;

	Grid(vector<std::tuple<double, double, double>> points, vector<Boundary> boundaries, 
				GridProperties properties, Eigen::VectorXd source, bool regularize);
	~Grid();

	void boundaryOp(std::string coarse);
	void setBCFlag(int boundary, std::string type, vector<double> boundValue);
	void build_laplacian();
	void push_inhomog_to_rhs();
	void build_normal_vecs(const char* filename, std::string geomtype);
	void fix_bounds_conn(std::vector<std::pair<int, int>>& ptsConn);
	void build_deriv_normal_bound();
	void sor(Eigen::SparseMatrix<double, 1>* matrix, Eigen::VectorXd* values, Eigen::VectorXd rhs);
	void setNeumannFlag();
	void modify_coeff_neumann(std::string coarse);
	void bound_eval_neumann();
	void fix_vector_bound_coarse (Eigen::VectorXd* vec);
	void print_bc_values();
	void print_bc_values(Eigen::VectorXd vec);
	void print_check_bc_normal_derivs();
	void build_interior_mat();

	double cond_L();
	Eigen::VectorXd residual();
	Eigen::VectorXd true_residual();

	void rcm_order_points();
	vector<Point> pointIDs_to_vector(const vector<int> &pointIDs);
	vector<int> kNearestNeighbors(Point point, bool neumannFlag, bool pointBCFlag, int k);
	vector <int> kNearestNeighbors(int pointNumber, bool neumannFlag, int k);

	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(Point point, bool neumann, bool pointBCFlag, int polyDeg);
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(int pointNum, bool neumann, int polyDeg);
	std::pair<Eigen::VectorXd, vector<int>> laplaceWeights(int pointID);
	std::pair<Eigen::VectorXd, vector<int>> derivx_weights(int pointID);
	std::pair<Eigen::VectorXd, vector<int>> derivy_weights(int pointID);
	std::pair<Eigen::VectorXd, vector<int>> pointInterpWeights(Point point, int polyDeg);
	
	int getSize();
	int getStencilSize();
	int getPolyDeg();
	vector<double> cond_rbf;
};
#endif
