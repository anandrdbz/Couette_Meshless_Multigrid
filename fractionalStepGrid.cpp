#include "fractionalStepGrid.hpp"
FractionalStepGrid::FractionalStepGrid(vector<std::tuple<double, double, double>> points, vector<Boundary> boundaries,
	GridProperties properties, Eigen::VectorXd source, bool regularize) : Grid(points, boundaries, properties, source, regularize){
	
	u = new Eigen::VectorXd(laplaceMatSize_);
	u_old = new Eigen::VectorXd(laplaceMatSize_);
	v = new Eigen::VectorXd(laplaceMatSize_);
	v_old = new Eigen::VectorXd(laplaceMatSize_);
	u_hat = new Eigen::VectorXd(laplaceMatSize_);
	v_hat = new Eigen::VectorXd(laplaceMatSize_);
	u_old->setZero();
	v_old->setZero();
	u->setZero();
	v->setZero();
	u_hat->setZero();
	v_hat->setZero();
}
FractionalStepGrid::~FractionalStepGrid() {
	delete u_hat;
	delete v_hat;
	delete u;
	delete v;
	delete u_old;
	delete v_old;
}
void FractionalStepGrid::prescribe_soln() {
	double x, y;
	double re = rho / mu;
	lambda = 0.5 * re - std::sqrt(0.25*re*re + 4*EIGEN_PI*EIGEN_PI);
	for (int i = 0; i < laplaceMatSize_; i++) {
		x = std::get<0>(points_[i]);
		y = std::get<1>(points_[i]);
		u->coeffRef(i) = 1 - std::exp(lambda * x) * std::cos(2 * EIGEN_PI * y);
		v->coeffRef(i) = lambda / (2 * EIGEN_PI) * std::exp(lambda * x) * std::sin(2 * EIGEN_PI * y);
		u_old->coeffRef(i) = 1 - std::exp(lambda * x) * std::cos(2 * EIGEN_PI * y);
		v_old->coeffRef(i) = lambda / (2 * EIGEN_PI) * std::exp(lambda * x) * std::sin(2 * EIGEN_PI * y);
		values_->coeffRef(i) = 0.5 * std::exp(2 * lambda * x);
	}
	values_->coeffRef(laplaceMatSize_) = 0;
}
void FractionalStepGrid::set_uv_bound() {
	double x, y,nx,ny;
	double re = rho / mu;
	if (flowType.compare("couette") == 0) {
		for (size_t i = 0; i < boundaries_.size(); i++) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				//cout << i << " " << j << " " << u->rows() << endl;
				x = std::get<0>(points_[boundaries_[i].bcPoints.at(j)]);
				y = std::get<1>(points_[boundaries_[i].bcPoints.at(j)]);
				nx = std::get<0>(normalVecs_[boundaries_[i].bcPoints.at(j)]);
				ny = std::get<1>(normalVecs_[boundaries_[i].bcPoints.at(j)]);
				if (std::abs(std::pow(x,2) + std::pow(y,2) - 1) < 1e-6 ) {
					u->coeffRef(boundaries_[i].bcPoints.at(j)) = 0;
					v->coeffRef(boundaries_[i].bcPoints.at(j)) = 0;
					u_old->coeffRef(boundaries_[i].bcPoints.at(j)) = 0;
					v_old->coeffRef(boundaries_[i].bcPoints.at(j)) = 0;
				}

				if (std::abs(std::pow(x,2) + std::pow(y,2) - 0.25) < 1e-6 ) {
					u->coeffRef(boundaries_[i].bcPoints.at(j)) = -ny;
					v->coeffRef(boundaries_[i].bcPoints.at(j)) = nx;
					u_old->coeffRef(boundaries_[i].bcPoints.at(j)) = -ny;
					v_old->coeffRef(boundaries_[i].bcPoints.at(j)) = nx;
				}


				/*if(j==0){
					cout<<v->coeffRef(boundaries_[0].bcPoints.at(0))<<endl;
					cout<<v->coeffRef(boundaries_[1].bcPoints.at(0))<<endl;
				}
				*/
			}
		}
		//cout<<v->coeffRef(boundaries_[0].bcPoints.at(0))<<endl;
		//cout<<v->coeffRef(boundaries_[1].bcPoints.at(1))<<endl;
	}
}
void FractionalStepGrid::build_derivX_mat() {
	vector<Eigen::Triplet<double>> tripletList;

	for (int i = 0; i < laplaceMatSize_; i++) {
		std::pair <Eigen::VectorXd, vector<int>> weights = derivx_weights(i);
		for (size_t j = 0; j < weights.second.size(); j++) {
			tripletList.push_back(Eigen::Triplet<double>(i, weights.second[j], weights.first(j))); 
		}
	}
	derivXMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(laplaceMatSize_, laplaceMatSize_);
	derivXMat_->setFromTriplets(tripletList.begin(), tripletList.end());
	derivXMat_->makeCompressed();

}

void FractionalStepGrid::build_derivY_mat() {
	vector<Eigen::Triplet<double>> tripletList;
	vector<Eigen::Triplet<double>> boundaryList;

	for (int i = 0; i < laplaceMatSize_; i++) {
		std::pair <Eigen::VectorXd, vector<int>> weights = derivy_weights(i);
		for (size_t j = 0; j < weights.second.size(); j++) {
			tripletList.push_back(Eigen::Triplet<double>(i, weights.second[j], weights.first(j)));
		}
	}
	derivYMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(laplaceMatSize_, laplaceMatSize_);
	derivYMat_->setFromTriplets(tripletList.begin(), tripletList.end());
	derivYMat_->makeCompressed();

	
}
void FractionalStepGrid::build_uv_laplace_mat() {
	vector<Eigen::Triplet<double>> tripletList;
	vector<Eigen::Triplet<double>> boundaryList;

	for (int i = 0; i < laplaceMatSize_; i++) {
		std::pair <Eigen::VectorXd, vector<int>> weights = laplaceWeights(i);
		for (size_t j = 0; j < weights.second.size(); j++) {
			tripletList.push_back(Eigen::Triplet<double>(i, weights.second[j], weights.first(j)));
		}
	}
	uvLaplaceMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(laplaceMatSize_, laplaceMatSize_);
	uvLaplaceMat_->setFromTriplets(tripletList.begin(), tripletList.end());
	uvLaplaceMat_->makeCompressed();

	
}
void FractionalStepGrid::calc_u_hat() {
	Eigen::VectorXd u_x, u_y, del2_u, p_x; 
	u_x = *derivXMat_ * *u;
	u_y = *derivYMat_ * *u;
	
	del2_u = *uvLaplaceMat_ * *u;

	
	for (int i = 0; i < laplaceMatSize_; i++) {
		u_hat->coeffRef(i) = u->coeff(i) + dt * (-(u->coeff(i) * u_x.coeff(i) + v->coeff(i) * u_y.coeff(i)) + mu / rho * del2_u.coeff(i));
	}


}
void FractionalStepGrid::calc_v_hat() {
	Eigen::VectorXd v_x, v_y, del2_v, p_y;
	v_x = *derivXMat_ * *v;
	v_y = *derivYMat_ * *v;

	del2_v = *uvLaplaceMat_ * *v;


	
	for (int i = 0; i < laplaceMatSize_; i++) {
		v_hat->coeffRef(i) = v->coeff(i) +dt * (-(u->coeff(i) * v_x.coeff(i) + v->coeff(i) * v_y.coeff(i)) + mu / rho * del2_v.coeff(i));
	}


}
void FractionalStepGrid::set_ppe_source() {
	//cout << laplaceMatSize_ << " " << source_.rows() << endl;
	source_.head(laplaceMatSize_) = rho/dt*(*derivXMat_ * *u_hat + *derivYMat_ * *v_hat); //interior
	//boundary
	double dpdx, dpdy, nx, ny;
	int currPoint;
	
	for (size_t i = 0; i < boundaries_.size(); i++) {
		for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
			currPoint = boundaries_[i].bcPoints[j];
			dpdx = -rho / dt * (u->coeff(currPoint) - u_hat->coeff(currPoint));
			dpdy = -rho / dt * (v->coeff(currPoint) - v_hat->coeff(currPoint));
			//cout << dpdx << " " << dpdy << endl;
			nx = std::get<0>(normalVecs_[currPoint]);
			ny = std::get<1>(normalVecs_[currPoint]);
			source_(currPoint) = nx * dpdx + ny * dpdy;
			//cout << source_(i) << endl;
		}
	}
	
}
void FractionalStepGrid::correct_u() { 
	*u = *u_hat - dt / rho * (*derivXMat_ * values_->head(laplaceMatSize_));
}
void FractionalStepGrid::correct_v() {
	*v = *v_hat - dt / rho * (*derivYMat_ * values_->head(laplaceMatSize_));
}
double FractionalStepGrid::fs_residual() {
	return (*u - *u_hat).lpNorm<1>()/laplaceMatSize_;
}
double FractionalStepGrid::delu() {
	return (*u - *u_old).lpNorm<1>()/laplaceMatSize_;
}