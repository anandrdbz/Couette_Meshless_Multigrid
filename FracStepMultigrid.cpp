#include "FracStepMultigrid.hpp"
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
FractionalStepMultigrid::FractionalStepMultigrid() {
	grids_ = vector<std::pair<int, FractionalStepGrid*>>(); //int = grid size
	residuals_ = vector<double>();
	residuals_.push_back(1.0);
}
FractionalStepMultigrid::~FractionalStepMultigrid() {
	for (size_t i = 0; i < grids_.size(); i++) {
		delete grids_.at(i).second;
		delete prolongMatrices_.at(i);
		delete restrictionMatrices_.at(i);
	}
}
Eigen::SparseMatrix<double>* FractionalStepMultigrid::buildInterpMatrix(FractionalStepGrid* baseGrid, FractionalStepGrid* targetGrid) {

	Eigen::SparseMatrix<double>* interpMatrix = new Eigen::SparseMatrix<double>(targetGrid->getSize(), baseGrid->getSize());
	vector<Eigen::Triplet<double>> tripletList;
	//tripletList: <row, col, value>
	std::pair<Eigen::VectorXd, vector<int>> pointWeights;
	for (int i = 0; i < targetGrid->getSize(); i++) {
		if(baseGrid->laplaceMatSize_ >= targetGrid->laplaceMatSize_){
		 	pointWeights = baseGrid->pointInterpWeights(targetGrid->points_[i], baseGrid->properties_.polyDeg);
		}
		else{
			pointWeights = baseGrid->pointInterpWeights(targetGrid->points_[i], targetGrid->properties_.polyDeg);
		}
		for (size_t j = 0; j < pointWeights.second.size(); j++) {

			tripletList.push_back(Eigen::Triplet<double>(i, pointWeights.second[j], pointWeights.first(j)));
		}
	}
	interpMatrix->setFromTriplets(tripletList.begin(), tripletList.end());
	return interpMatrix;
}

void FractionalStepMultigrid::buildProlongMatrices() {
	prolongMatrices_.resize(grids_.size());
	for (size_t i = 0; i < grids_.size() - 1; i++) {
		prolongMatrices_[i] = (buildInterpMatrix(grids_[i].second, grids_[i + 1].second));
	}
	prolongMatrices_[prolongMatrices_.size() - 1] = NULL;
}
void FractionalStepMultigrid::buildRestrictionMatrices() {
	restrictionMatrices_.resize(grids_.size());
	restrictionMatrices_[0] = NULL;
	for (size_t i = 1; i < grids_.size(); i++) {
		restrictionMatrices_[i] = buildInterpMatrix(grids_[i].second, grids_[i - 1].second);
	}
}
void FractionalStepMultigrid::buildMatrices() {
	buildProlongMatrices();
	buildRestrictionMatrices();
	Grid* currGrid;
	//set coarse source coeffs back to 0
	for (size_t i = 0; i < grids_.size() - 1; i++) {
		currGrid = grids_[i].second;
		if (i != grids_.size() - 1) {
			currGrid->modify_coeff_neumann("coarse");
		}
	}
}

void FractionalStepMultigrid::vCycle() {
	//Restriction
	Grid* currGrid;
	currGrid = grids_[grids_.size() - 1].second;
	if (grids_.size() == 1) {
		currGrid->sor(currGrid->laplaceMat_, currGrid->values_, (currGrid->source_));
		return;
	}

	currGrid->bound_eval_neumann();
	//std::cout << std::setprecision(12) << "PPE Residual: " << resid_norm << std::endl;
	//Restriction
	for (size_t i = grids_.size() - 1; i > 0; i--) {
		currGrid = grids_[i].second;
		std::string gridType = i == grids_.size() - 1 ? "fine" : "coarse";
		//cout << gridType << endl;
		if (i != grids_.size() - 1) {
			currGrid->values_->setZero();
		}
		currGrid->boundaryOp(gridType);
		currGrid->sor(currGrid->laplaceMat_, currGrid->values_, (currGrid->source_));

		(grids_[i - 1].second->source_)(Eigen::seq(0, grids_[i - 1].second->laplaceMatSize_ - 1)) = (*(restrictionMatrices_[i])) * (currGrid->residual())(Eigen::seq(0, currGrid->laplaceMatSize_ - 1));
		grids_[i - 1].second->fix_vector_bound_coarse(&grids_[i - 1].second->source_);
		if (currGrid->regularizeFlag_) {
			grids_[i - 1].second->source_.coeffRef(grids_[i - 1].second->source_.rows() - 1) = 0;
		}

	}
	//cout << grids_[0].second->source_.norm() << endl;
	//Iterate on coarsest grid
	currGrid->boundaryOp("coarse");
	currGrid = grids_[0].second;
	currGrid->values_->setZero();
	currGrid->sor(currGrid->laplaceMat_, currGrid->values_, (currGrid->source_));
	currGrid->sor(currGrid->laplaceMat_, currGrid->values_, (currGrid->source_));

	//correction and prolongation
	Eigen::VectorXd correction;
	for (size_t i = 1; i < grids_.size(); i++) {
		currGrid = grids_[i].second;
		//correction
		correction = (*prolongMatrices_[i - 1]) * ((*grids_[i - 1].second->values_)(Eigen::seq(0, grids_[i - 1].second->laplaceMatSize_ - 1)));
		currGrid->fix_vector_bound_coarse(&correction);

		(*currGrid->values_)(Eigen::seq(0, currGrid->laplaceMatSize_ - 1)) += correction;
		//smoother
		currGrid->sor(currGrid->laplaceMat_, currGrid->values_, (currGrid->source_));
	}
	currGrid->boundaryOp("fine");
}

double FractionalStepMultigrid::gmres(double tol, std::string Solver){
	Grid* finestGrid;
	finestGrid = grids_[grids_.size()-1].second;

	Eigen::VectorXd r = finestGrid->rhs_ - *(finestGrid->interior_mat) * (*finestGrid->sol_);
	Eigen::VectorXd z,p,w,v;


	double alpha,l2res,l2res_prev;
	l2res = 0;
	int kmax;
	vector<Eigen::VectorXd> P;
	vector<Eigen::VectorXd> W;

	kmax = finestGrid->laplaceMatSize_;

	
	for(int k = 0; k < kmax; k++){

		//GMRES with Multigrid
		if(Solver == "Multigrid"){
			finestGrid->source_ = *(finestGrid->prolong_mat) * r;
			finestGrid->values_->setZero();
			for(int j = 0; j < 3; j++){
				vCycle();
				//double mgtol = finestGrid->residual().lpNorm<1>() / finestGrid->source_.lpNorm<1>(); 
				//cout << "mgtol " << mgtol << endl;

			}
			
			z = *(finestGrid->restrict_mat) * *(finestGrid->values_);			
		}
		else if(Solver == "ILU"){
			z = finestGrid->solver.solve(r);
		}
		else if(Solver == "SOR"){
			finestGrid->source_ = *(finestGrid->prolong_mat) * r;
			finestGrid->values_->setZero();
			finestGrid->sor(finestGrid->laplaceMat_,finestGrid->values_,(finestGrid->source_));
			//finestGrid->sor(finestGrid->laplaceMat_,finestGrid->values_,&(finestGrid->source_));
			z = *(finestGrid->restrict_mat) * *(finestGrid->values_);
		}
        else if(Solver == "None"){
			z = r;        	
        }

		if(k == 0){
			p = z;
			P.push_back(p);
			w = *(finestGrid->interior_mat) * p;
			W.push_back(w);
		}
		else{
			p = z;
			v = *(finestGrid->interior_mat) * z;
			for(int i = 0; i<k; i++){
				double beta = W[i].dot(v)/(W[i].dot(W[i]));
				p = p - beta*P[i];
			}
			P.push_back(p);
			w = *(finestGrid->interior_mat) * p;
			W.push_back(w);
		}
		alpha = w.dot(r)/(w.dot(w));
		*(finestGrid->sol_) = *(finestGrid->sol_) + alpha*p;
		r = r - alpha*w;
		Eigen::VectorXd b = finestGrid->rhs_;
		l2res_prev = l2res;
		l2res = r.lpNorm<2>() / (b.lpNorm<2>()) ;

		residuals_.push_back(l2res);
		if(l2res < tol){
			cout << "l2res: " << l2res << endl;
			break;
		}
		/*if(std::abs(l2res_prev - l2res) < 1e-12){
			break;
		}*/
	}
	return l2res;

}



void FractionalStepMultigrid::solveLoop() {

}
double FractionalStepMultigrid::residual() {
	Grid* finegrid = grids_[grids_.size() - 1].second;
	return finegrid->residual().lpNorm<1>() / finegrid->source_.lpNorm<1>();
}
void FractionalStepMultigrid::addGrid(FractionalStepGrid* grid) {
	grids_.push_back(std::pair<int, FractionalStepGrid*>(grid->getSize(), grid));
	sortGridsBySize();
}
void FractionalStepMultigrid::sortGridsBySize() {
	std::sort(grids_.begin(), grids_.end());
}
