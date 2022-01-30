#ifndef GRID_CLASSES_H
#define GRID_CLASSES_H
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Eigen/Core"
class GridProperties {
public:
	int rbfExp;
	int polyDeg;
	int laplaceMatSize;
	int stencilSize;
	double omega;
	int iters;
};
class Boundary {
public:
	int type;
	std::vector<int> bcPoints;
	std::vector<double> values;
};
class deriv_normal_bc {
public:
	int pointID;
	Eigen::VectorXd weights;
	std::vector<int> neighbors;
	double value;

};
#endif