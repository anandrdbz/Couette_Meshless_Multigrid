#ifndef MULTIGRID_H
#define MULTIGRID_H
#include "grid.h"
#include <Eigen/Core>
class Multigrid {
public:
	vector<std::pair<int, Grid*>> grids_;
	vector<int> sorGridIters_;
	vector<Eigen::SparseMatrix<double>*> restrictionMatrices_;
	vector<Eigen::SparseMatrix<double>*> prolongMatrices_;
	vector<double> residuals_;

	void sortGridsBySize();
	Eigen::SparseMatrix<double>* buildInterpMatrix(Grid* baseGrid, Grid* targetGrid);
	void buildRestrictionMatrices();
	void buildProlongMatrices();

	Multigrid();
	~Multigrid();
	void addGrid(Grid* grid);
	void buildMatrices();
	void vCycle();
	double residual();
};
#endif