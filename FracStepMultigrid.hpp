#ifndef FRAC_STEP_MG_H
#define FRAC_STEP_MG_H
#include "fractionalStepGrid.hpp"
class FractionalStepMultigrid{
public:
	vector<std::pair<int, FractionalStepGrid*>> grids_;
	vector<int> sorGridIters_;
	vector<Eigen::SparseMatrix<double>*> restrictionMatrices_;
	vector<Eigen::SparseMatrix<double>*> prolongMatrices_;
	vector<double> residuals_;

	void sortGridsBySize();
	Eigen::SparseMatrix<double>* buildInterpMatrix(FractionalStepGrid* baseGrid, FractionalStepGrid* targetGrid);
	void buildRestrictionMatrices();
	void buildProlongMatrices();

	FractionalStepMultigrid();
	~FractionalStepMultigrid();
	void addGrid(FractionalStepGrid* grid);
	void buildMatrices();
	double gmres(double tol, std::string Solver);
	void vCycle();
	void solveLoop();
	double residual();

};
#endif