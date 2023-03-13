#include "grid.h"     
#define pi 3.141592653589793238462643383279

using std::cout;
using std::endl;
Grid::Grid(vector<Point> points, vector<Boundary> boundaries,
	GridProperties properties, Eigen::VectorXd source, bool regularize) {
	int numPoint = (int)(points.size());
	points_ = points;
	boundaries_ = boundaries;
	properties_ = properties;
	source_ = source;
	cond_rbf = vector<double>();
	neumannFlag_ = false;
	regularizeFlag_ = regularize;
	setNeumannFlag();

	laplaceMatSize_ = numPoint;
	int A_size = numPoint + (int) regularizeFlag_;
	bcFlags_ = std::vector<int>(numPoint + (int) regularizeFlag_, 0);
	normalVecs_ = std::vector<Point>(numPoint);
	laplaceMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(A_size, A_size);
	laplaceMat_->setZero();
	values_ = new Eigen::VectorXd(A_size);
	values_->setZero();
	
	bcount_ = 0;
	for (size_t i = 0; i < boundaries_.size(); i++) {
		bcount_ += (boundaries_[i].bcPoints).size();
	}

	sol_ = new Eigen::VectorXd(A_size - bcount_);
	sol_->setZero();
	rhs_ = Eigen::VectorXd(A_size - bcount_);
	rhs_.setZero();
	
	restrict_mat = new Eigen::SparseMatrix<double, Eigen::RowMajor>(A_size - bcount_, A_size);
	restrict_mat->setZero();
	

	prolong_mat = new Eigen::SparseMatrix<double, Eigen::RowMajor>(A_size , A_size - bcount_);
	prolong_mat->setZero();

	interior_mat = new Eigen::SparseMatrix<double, Eigen::RowMajor>(A_size - bcount_, A_size - bcount_);
	interior_mat->setZero();

	neumann_boundary_coeffs_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(A_size, A_size);
	diags = Eigen::VectorXd(A_size);
}
Grid::~Grid() {
	delete values_;
	delete laplaceMat_;
}

void Grid::setBCFlag(int bNum, std::string type, vector<double> boundValues) {
	for(int i = 0; i < laplaceMatSize_+(int) regularizeFlag_; i++){
		bcFlags_[i] = 0;
	}
	for (size_t i = 0; i < boundaries_.size(); i++) {
		for(size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++){
			bcFlags_[(boundaries_[i].bcPoints).at(j)] = boundaries_[i].type;
		}
	}
}

void Grid::boundaryOp(std::string coarse) {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		//dirichlet
		if ((boundaries_[i]).type == 1) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				(*values_)((boundaries_[i].bcPoints).at(j)) = coarse.compare("coarse") == 0 ? 0 : (boundaries_[i].values).at(j);
			}
		}
	}
}
void Grid::setNeumannFlag() {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		if (boundaries_[i].type == 2) {
			neumannFlag_ = true;
			return;
		}
	}
	neumannFlag_ = false;
}

void Grid::modify_coeff_neumann(std::string coarse) {

	for (size_t i = 0; i < boundaries_.size(); i++) {
		if ((boundaries_[i]).type == 2) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				source_(boundaries_[i].bcPoints.at(j)) = coarse.compare("coarse") == 0 ? 0 : boundaries_[i].values.at(j);
			}
		}
	}
	if(regularizeFlag_){
		source_.coeffRef(source_.rows() - 1) = 0;
	}
	
}
void Grid::bound_eval_neumann() {
	const double* laplaceValues = laplaceMat_->valuePtr();
	//column indices of values
	const int* innerValues = laplaceMat_->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = laplaceMat_->outerIndexPtr();
	int rowStart, rowEnd, curr;
	double diag, boundValue;
	for (size_t i = 0; i < boundaries_.size(); i++) {
		if ((boundaries_[i]).type == 2) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				curr = boundaries_[i].bcPoints[j];
				//values_->coeffRef(curr) = std::cos(pi*std::get<0>(points_[curr]))*std::cos(pi*std::get<1>(points_[curr]));

				rowStart = outerValues[curr];
				rowEnd = outerValues[curr + 1];
				boundValue = source_.coeff(curr);
				for (int k = rowStart; k < rowEnd; k++) {
					if (innerValues[k] == curr) {
						diag = laplaceValues[k];
						continue;
					}
					boundValue -= values_->coeff(innerValues[k]) * laplaceValues[k];
				}
				boundValue /= diag;
				values_->coeffRef(curr) = boundValue;

			}
		}
	}
}
void Grid::sor(Eigen::SparseMatrix<double, Eigen::RowMajor>* matrix, Eigen::VectorXd* values, Eigen::VectorXd rhs) {
	double x_i, diagCoeff;
	const double* laplaceValues = matrix->valuePtr();
	const int* innerValues = matrix->innerIndexPtr();
	const int* outerValues = matrix->outerIndexPtr();
	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	int numNonZeros = matrix->nonZeros();
	int inner;
	for (int it = 0; it < properties_.iters; it++) {
		valueIdx = 0;
		innerIdx = 0;
		outerIdx = 0;

		for (int i = 0; i < matrix->rows(); i++) {
			if (!(regularizeFlag_ && i == matrix->rows() - 1) && bcFlags_[i] != 0) {
				outerIdx++;
				continue;
			}
			x_i = 0;
			diagCoeff = 0;
			rowStartIdx = outerValues[outerIdx];
			rowEndIdx = outerValues[outerIdx + 1];
			for (int j = rowStartIdx; j < rowEndIdx; j++) {

				if (innerValues[j] == i) {
					diagCoeff = laplaceValues[j];
					continue;
				}


				inner = innerValues[j];
				x_i -= laplaceValues[j] * values->coeff(inner);
			}
			x_i += rhs.coeff(i);

			x_i *= properties_.omega / diagCoeff;
			x_i += (1 - properties_.omega) * values->coeff(i);
			values->coeffRef(i) = x_i;
			outerIdx++;
		}
		bound_eval_neumann();
	}
}
Eigen::VectorXd Grid::residual() {
	Eigen::VectorXd res = source_ - (*laplaceMat_) * (*values_);
	fix_vector_bound_coarse(&res);
	return res;
}

Eigen::VectorXd Grid::true_residual(){
	Eigen::VectorXd res = rhs_ - (*interior_mat) * (*sol_);
	return res;
}

double Grid::cond_L() {
	double cond = 1.0 / laplaceMat_->toDense().partialPivLu().rcond();
	return cond;
}
void Grid::print_bc_values(Eigen::VectorXd vec) {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		if ((boundaries_[i]).type == 1) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				cout << "bc value: " << vec.coeff(boundaries_[i].bcPoints.at(j)) << endl;
			}
		}
	}
}
void Grid::print_bc_values() {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
			cout << "bc value: " << boundaries_[i].values.at(j) << endl;
		}
	}
}
void Grid::print_check_bc_normal_derivs() {
	double* laplaceValues = laplaceMat_->valuePtr();
	const int* innerValues = laplaceMat_->innerIndexPtr();
	const int* outerValues = laplaceMat_->outerIndexPtr();
	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	double sum;
	valueIdx = 0;
	innerIdx = 0;
	outerIdx = 0;

	for (int i = 0; i < laplaceMat_->rows() - (int) regularizeFlag_; i++) {
		if (bcFlags_[i] != 2) {
			outerIdx++;
			continue;
		}
		sum = 0;
		rowStartIdx = outerValues[outerIdx];
		rowEndIdx = outerValues[outerIdx + 1];
		for (int j = rowStartIdx; j < rowEndIdx; j++) {
			sum += laplaceValues[j] * values_->coeff(innerValues[j]);
		}
		cout << "d/dn: " << sum << " rhs: " << source_.coeff(i) << endl;
		outerIdx++;
	}
}
void Grid::fix_vector_bound_coarse(Eigen::VectorXd* vec) {
	for (size_t i = 0; i < boundaries_.size(); i++) {

		for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
			vec->coeffRef(boundaries_[i].bcPoints.at(j)) = 0;
		}
		
	}
}
vector<Point> Grid::pointIDs_to_vector(const vector<int>& pointIDs) {
	vector<Point> points;
	for (size_t i = 0; i < pointIDs.size(); i++) {
		points.push_back(points_[pointIDs.at(i)]);
	}
	return points;
}
vector<int> Grid::kNearestNeighbors(int pointID, bool neumann, int stencilSize) {
	return kNearestNeighbors(points_[pointID], neumann, (bcFlags_[pointID] != 0), stencilSize);
}
vector<int> Grid::kNearestNeighbors(Point refPoint, bool neumann, bool pointBCFlag, int stencilSize) {
	vector<std::pair<double, int>> distances;
	Point queryPoint;
	int samePoint = -1;
	for (int i = 0; i < laplaceMatSize_; i++) {
		queryPoint = points_[i];
		distances.push_back(std::pair<double, int>(distance(refPoint, queryPoint), i));

		if (distances[i].first == 0) {
			samePoint = i;
		}

	}

	//max-heap selection algorithm
	int k = stencilSize;

	int lastInitIndex = k;
	vector<std::pair<double, int>> maxHeap;
	for (int i = 0; i < lastInitIndex; i++) {
		if (i != samePoint && (pointBCFlag && neumann && bcFlags_[i] != 0)) {
			lastInitIndex++;
			continue;
		}
		maxHeap.push_back(distances[i]);
	}
	std::make_heap(maxHeap.begin(), maxHeap.end());
	for (int i = lastInitIndex; i < laplaceMatSize_; i++) {
		if (i == samePoint || !(pointBCFlag && neumann && bcFlags_[i] != 0) && distances[i] < maxHeap.front()) {
			std::pop_heap(maxHeap.begin(), maxHeap.end());
			maxHeap.pop_back();

			maxHeap.push_back(distances[i]);
			std::push_heap(maxHeap.begin(), maxHeap.end());
		}
	}
	std::sort_heap(maxHeap.begin(), maxHeap.end());
	vector<int> nearestNeighbors;
	//Have to include point itself in the stencil since otherwise diag would be zeros.
	for (int i = 0; i < k; i++) {
		nearestNeighbors.push_back(maxHeap[i].second);
	}
	return nearestNeighbors;

}
//Builds PHS coefficients matrix for a point
//returns pair with first = coeff matrix and second = list of neighbors, to save time.
std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> Grid::buildCoeffMatrix(Point point, bool neumann, bool pointBCFlag, int polyDeg) {
	//cout << properties_.polyDeg << endl;
	//cout << polyDeg << endl;
	int polyTerms = (polyDeg + 1) * (polyDeg + 2) / 2;
	int stencilSize = (int)(2.0 * (polyDeg + 1) * (polyDeg + 2) / 2);
	vector<int> neighbors = kNearestNeighbors(point, neumann, pointBCFlag, stencilSize);
	vector<Point> scaledPoints = shifting_scaling(pointIDs_to_vector(neighbors), point);
	//build A-matrix
	Eigen::MatrixXd coeff_mat = Eigen::MatrixXd::Zero(stencilSize + polyTerms, stencilSize + polyTerms);

	double a_coeff, r;
	for (int i = 0; i < stencilSize; i++) {
		for (int j = i; j < stencilSize; j++) {
			r = distance(scaledPoints[i], scaledPoints[j]);
			a_coeff = std::pow(r, properties_.rbfExp);
			coeff_mat(i, j) = a_coeff;
			coeff_mat(j, i) = a_coeff;
		}
	}
	//fill P-matrix sections.
	int colIndex;
	double x, y, p_coeff;
	for (int row = 0; row < stencilSize; row++) {
		colIndex = stencilSize;
		for (int p = 0; p <= polyDeg; p++) {
			for (int q = 0; q <= p; q++) {
				x = std::get<0>(scaledPoints[row]);
				y = std::get<1>(scaledPoints[row]);
				p_coeff = std::pow(x, p - q) * std::pow(y, q);
				coeff_mat(row, colIndex) = p_coeff;
				coeff_mat(colIndex, row) = p_coeff;
				colIndex++;
			}
		}
	}
	return std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>>(coeff_mat, neighbors, scaledPoints);
}

std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> Grid::buildCoeffMatrix(int pointID, bool neumann, int polyDeg) {
	return buildCoeffMatrix(points_[pointID], neumann, (bcFlags_[pointID] != 0), polyDeg);
}
std::pair<Eigen::VectorXd, vector<int>> Grid::derivx_weights(int pointID) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(pointID, neumannFlag_, properties_.polyDeg);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	int polyTerms = (properties_.polyDeg + 1) * (properties_.polyDeg + 2) / 2;

	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);

	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	double yRef, xRef;
	double M = (double)(properties_.rbfExp);
	double coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		xRef = std::get<0>(scaledPoints[i]);
		yRef = std::get<1>(scaledPoints[i]);
		if (i > 0) {
			rhs(i) = M * std::pow(distance(scaledPoints[i], evalPoint), M - 2) * (xEval - xRef);
		}
	}
	int rowIndex = properties_.stencilSize;
	for (int p = 0; p <= properties_.polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			double laplacePoly = 0;
			if (p - q - 1 >= 0) {
				laplacePoly += (p - q) * std::pow(xEval, p - q - 1) * std::pow(yEval, q);
			}
			rhs(rowIndex) = laplacePoly;
			rowIndex++;
		}
	}
	Eigen::VectorXd weights = std::get<0>(coeffs).fullPivLu().solve(rhs);

	double scale = std::get<0>(scaledPoints[scaledPoints.size() - 2]);
	for (int i = 0; i < weights.rows(); i++) {
		weights(i) /= scale;
	}
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
std::pair<Eigen::VectorXd, vector<int>> Grid::derivy_weights(int pointID) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(pointID, neumannFlag_, properties_.polyDeg);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	int polyTerms = (properties_.polyDeg + 1) * (properties_.polyDeg + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);

	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	double yRef, xRef;
	double M = (double)(properties_.rbfExp);
	double coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		xRef = std::get<0>(scaledPoints[i]);
		yRef = std::get<1>(scaledPoints[i]);
		if (i > 0) {
			rhs(i) = M * std::pow(distance(scaledPoints[i], evalPoint), M - 2) * (yEval - yRef);
		}
	}
	int rowIndex = properties_.stencilSize;
	for (int p = 0; p <= properties_.polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			double laplacePoly = 0;
			if (q - 1 >= 0) {
				laplacePoly += q * std::pow(xEval, p - q) * std::pow(yEval, q - 1);
			}
			rhs(rowIndex) = laplacePoly;
			rowIndex++;
		}
	}
	Eigen::VectorXd weights = std::get<0>(coeffs).fullPivLu().solve(rhs);
	double scale = std::get<0>(scaledPoints[scaledPoints.size() - 2]);
	for (int i = 0; i < weights.rows(); i++) {
		weights(i) /= scale;
	}
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
std::pair<Eigen::VectorXd, vector<int>> Grid::laplaceWeights(int pointID) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(pointID, neumannFlag_, properties_.polyDeg);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	int polyTerms = (properties_.polyDeg + 1) * (properties_.polyDeg + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);

	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	double yRef, xRef, D;
	double M = (double)(properties_.rbfExp);
	double coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		xRef = std::get<0>(scaledPoints[i]);
		yRef = std::get<1>(scaledPoints[i]);
		D = (xEval * xEval - 2 * xEval * xRef + xRef * xRef + yEval * yEval - 2 * yEval * yRef + yRef * yRef);
		if (D > 0) {
			rhs(i) = (std::pow(2 * xEval - 2 * xRef, 2) + std::pow(2 * yEval - 2 * yRef, 2)) * (M / 2) * (M / 2 - 1) * std::pow(D, M / 2 - 2)
				+ 2 * M * std::pow(D, M / 2 - 1);
		}
	}

	int rowIndex = properties_.stencilSize;
	for (int p = 0; p <= properties_.polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			double laplacePoly = 0;
			if (p - q - 2 >= 0) {
				laplacePoly += (p - q) * (p - q - 1) * std::pow(xEval, p - q - 2) * std::pow(yEval, q);
			}
			if (q - 2 >= 0) {
				laplacePoly += q * (q - 1) * std::pow(xEval, p - q) * std::pow(yEval, q - 2);
			}
			rhs(rowIndex) = laplacePoly;
			rowIndex++;
		}
	}
	Eigen::VectorXd weights = std::get<0>(coeffs).fullPivLu().solve(rhs);
	double scale = std::get<0>(scaledPoints[scaledPoints.size() - 2]);
	for (int i = 0; i < weights.rows(); i++) {
		weights(i) /= std::pow(scale, 2);
	}
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
void Grid::fix_bounds_conn(std::vector<std::pair<int, int>>& ptsConn) {
	int start_bc_point, curr_bc_point, next_bc_point, temp;
	for (size_t b = 0; b < boundaries_.size(); b++) {
		start_bc_point = boundaries_[b].bcPoints[0];
		curr_bc_point = start_bc_point;
		do {
			next_bc_point = ptsConn[curr_bc_point].first;

			if (ptsConn[next_bc_point].second != curr_bc_point) {
				temp = ptsConn[next_bc_point].first;
				ptsConn[next_bc_point].first = ptsConn[next_bc_point].second;
				ptsConn[next_bc_point].second = temp;
			}
			curr_bc_point = next_bc_point;
		} while (curr_bc_point != start_bc_point);
	}
}
void Grid::build_normal_vecs(const char* filename, std::string geomtype) {
	Point curr;
	double x, y,x0,x1,normPoint,nx,ny;



	for (size_t b = 0; b < boundaries_[0].bcPoints.size(); b++) {
		curr = points_[boundaries_[0].bcPoints[b]];
		x = std::get<0>(curr);
		y = std::get<1>(curr);
		normPoint = std::sqrt(x * x + y * y);
		nx = x/normPoint;
		ny = y/normPoint;
		normalVecs_[boundaries_[0].bcPoints[b]] = std::make_tuple(nx, ny, 0);
		

	}
/*
		for (size_t b = 0; b < boundaries_[1].bcPoints.size(); b++) {
		curr = points_[boundaries_[1].bcPoints[b]];
		x = std::get<0>(curr);
		y = std::get<1>(curr);
		x = x ;
		y = y ;
		normPoint = std::sqrt(x * x + y * y);
		x /= normPoint;
		y /= normPoint;
		normalVecs_[boundaries_[1].bcPoints[b]] = std::make_tuple(x, y, 0);
	}
*/					

}

void Grid::build_deriv_normal_bound() {
	double x, y, xWeight, yWeight;
	std::pair<Eigen::VectorXd, vector<int>> coeffs;
	int currPoint;
	deriv_normal_coeffs_ = vector<deriv_normal_bc>();
	deriv_normal_bc bound;
	double root2 = std::sqrt(2);
	for (size_t i = 0; i < boundaries_.size(); i++) {
		if ((boundaries_[i]).type == 2) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				currPoint = boundaries_[i].bcPoints[j];
				x = std::get<0>(points_[currPoint]);
				y = std::get<1>(points_[currPoint]);
				xWeight = std::get<0>(normalVecs_[currPoint]);
				yWeight = std::get<1>(normalVecs_[currPoint]);
				coeffs = derivx_weights(currPoint);
				coeffs.first *= xWeight;
				coeffs.first += yWeight * derivy_weights(currPoint).first;

				bound.pointID = currPoint;
				bound.value = boundaries_[i].values[j];
				bound.weights = coeffs.first;
				bound.neighbors = coeffs.second;
				deriv_normal_coeffs_.push_back(bound);
			}
		}
	}

}
void Grid::build_laplacian() {
	vector<Eigen::Triplet<double>> tripletList;
	vector<Eigen::Triplet<double>> boundaryList;

	for (int i = 0; i < laplaceMatSize_; i++) {
		std::pair <Eigen::VectorXd, vector<int>> weights = laplaceWeights(i);
		if (bcFlags_[i] != 2) {
			for (size_t j = 0; j < weights.second.size(); j++) {
				tripletList.push_back(Eigen::Triplet<double>(i, weights.second[j], weights.first(j)));
				if (bcFlags_[i] == 0 && bcFlags_[weights.second[j]] == 2) {
					boundaryList.push_back(Eigen::Triplet<double>(i, weights.second[j], weights.first(j)));
				}
				if (i == weights.second[j]) {
					diags.coeffRef(i) = weights.first(j);
				}
			}
		}
		if (regularizeFlag_ && bcFlags_[i] != 2) {
			tripletList.push_back(Eigen::Triplet<double>(i, laplaceMatSize_, 1));
		}
	}
	if (regularizeFlag_) {
		for (int i = 0; i < laplaceMatSize_ + (int) regularizeFlag_; i++) {
			if ((i == laplaceMatSize_) || bcFlags_[i] != 2) {
				tripletList.push_back(Eigen::Triplet<double>(laplaceMatSize_, i, 1));
			}
		}
	}
	//Add neumann boundary rows: d/dn = rhs
	deriv_normal_bc bound;
	if (neumannFlag_) {
		for (size_t i = 0; i < deriv_normal_coeffs_.size(); i++) {
			bound = deriv_normal_coeffs_[i];
			for (size_t j = 0; j < bound.neighbors.size(); j++) {
				tripletList.push_back(Eigen::Triplet<double>(bound.pointID, bound.neighbors[j], bound.weights.coeff(j)));
				if (bound.pointID == bound.neighbors[j]) {
					diags.coeffRef(bound.pointID) = bound.weights.coeff(j);
				}
			}
		}
	}
	laplaceMat_->setFromTriplets(tripletList.begin(), tripletList.end());
	laplaceMat_->makeCompressed();
	neumann_boundary_coeffs_->setFromTriplets(boundaryList.begin(), boundaryList.end());
	neumann_boundary_coeffs_->makeCompressed();
	if (!implicitFlag_) {
		return;
	}

	double* laplaceValues = laplaceMat_->valuePtr();
	const int* innerValues = laplaceMat_->innerIndexPtr();
	const int* outerValues = laplaceMat_->outerIndexPtr();
	int outerIdx_int, rowStartIdx_int, rowEndIdx_int,
		outerIdx_bound, rowStartIdx_bound, rowEndIdx_bound, j_col;
	double A_ij, A_jj, A_ik, A_jk, A_ij_mod;
	outerIdx_int = 0;
	int start,end;
	int ibflag;
	vector<int> iblist;


	for(int i = 0; i < laplaceMatSize_; i++){
		if (bcFlags_[i] != 0) {
			continue;
		}
		start = outerValues[i];
		end = outerValues[i+1];
		for(int j = start; j < end; j++){
			if(bcFlags_[innerValues[j]] != 0) {
				ibflag = 1; 
				break;
			}
		}		
		if(ibflag == 1){ 
			iblist.push_back(i);
			ibflag = 0;
		}
	}



	vector<std::pair<int, double>> tempRowValuesInt, tempRowValuesBound;
	for (int q = 0; q < iblist.size(); q++) {
		int i = iblist[q];
		rowStartIdx_int = outerValues[i];
		rowEndIdx_int = outerValues[i + 1];
		tempRowValuesInt.clear();
		//Iterate over row I to find all boundary points J.
		for (int j = rowStartIdx_int; j < rowEndIdx_int; j++) {
			if (innerValues[j] != laplaceMatSize_ && bcFlags_[innerValues[j]] != 0) {
				tempRowValuesInt.push_back(std::make_pair(innerValues[j], laplaceValues[j]));
			}
		}

		//iterate over every boundary point J.
		for (int j = 0; j < tempRowValuesInt.size(); j++) {
			//want to now add all stencil points in the row J.
			j_col = tempRowValuesInt[j].first;
			outerIdx_bound = j_col;
			A_ij = tempRowValuesInt[j].second;
			A_ij_mod = A_ij;

			rowStartIdx_bound = outerValues[outerIdx_bound];
			rowEndIdx_bound = outerValues[outerIdx_bound + 1];
			//find diagonal
			/*
			for (int k = rowStartIdx_bound; k < rowEndIdx_bound; k++) {
				if (innerValues[k] == j_col) {
					A_jj = laplaceValues[k];
					break;
				}
			}
			*/
			A_jj = diags.coeff(j_col);
			//Iterate over all of the stencil of boundary point J, A_jk.
			//We are iterating over the correct points, I did bcflag check.
			for (int k = rowStartIdx_bound; k < rowEndIdx_bound; k++) {
				A_jk = laplaceValues[k];
				if (j_col == innerValues[k]) {
					//A_ij_mod -= A_jk * A_ij / A_jj;
					continue;
				}
				tripletList.push_back(Eigen::Triplet<double>(i, innerValues[k], -A_jk * A_ij / A_jj));
			}
			//Decouple Boundary
			tripletList.push_back(Eigen::Triplet<double>(i, j_col, -A_ij));
			//cout << "bcI: " << bcFlags_[i] << "bcJ: " << bcFlags_[j_col] << "coeff: " << laplaceMat_->coeff(i, j_col) - A_ij << endl;
		}
	}
	delete laplaceMat_;
	laplaceMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(points_.size() + (int) regularizeFlag_, points_.size() + (int) regularizeFlag_);
	laplaceMat_->setFromTriplets(tripletList.begin(), tripletList.end());
	laplaceMat_->makeCompressed();

}
void Grid::push_inhomog_to_rhs() {
	if (!implicitFlag_) {
		return;
	}
	double* laplaceValues_bound = neumann_boundary_coeffs_->valuePtr();
	const int* innerValues = neumann_boundary_coeffs_->innerIndexPtr();
	const int* outerValues = neumann_boundary_coeffs_->outerIndexPtr();
	Eigen::VectorXd sourceCopy = source_;
	double diag, A_ij;
	for (int i = 0; i < laplaceMatSize_; i++) {
		if (bcFlags_[i] != 0) {
			continue;
		}
		//cout << i << endl;
		for (int j = outerValues[i]; j < outerValues[i + 1]; j++) {
			diag = diags.coeff(innerValues[j]);
			A_ij = laplaceValues_bound[j];
			source_(i) -= A_ij * sourceCopy.coeff(innerValues[j]) / diag;
		}
	}

}

std::pair<Eigen::VectorXd, vector<int>> Grid::pointInterpWeights(Point point, int polyDeg) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(point, false, false, polyDeg);
	int polyTerms = (polyDeg + 1) * (polyDeg + 2) / 2;
	int stencilSize = (int)(2.0 * (polyDeg + 1) * (polyDeg + 2) / 2);
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(stencilSize + polyTerms);
	Eigen::MatrixXd coeff_mat = std::get<0>(coeffs);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	//Build RHS rbf terms	
	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	for (int i = 0; i < stencilSize; i++) {
		rhs(i) = std::pow(distance(evalPoint, scaledPoints[i]), properties_.rbfExp);
	}
	//RHS poly terms
	int rowIndex = stencilSize;
	for (int p = 0; p <= polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			rhs(rowIndex) = std::pow(xEval, p - q) * std::pow(yEval, q);
			rowIndex++;
		}
	}
	Eigen::VectorXd weights = coeff_mat.fullPivLu().solve(rhs);
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
void Grid::rcm_order_points() {
	vector<vector<int>> adjacency;
	vector<int> order(points_.size());

	vector<int> neighbors;
	for (size_t i = 0; i < points_.size(); i++) {
		neighbors = kNearestNeighbors(points_[i], neumannFlag_, (bcFlags_[i] != 0), properties_.stencilSize);
		adjacency.push_back(neighbors);
	}
	if (neumannFlag_ && implicitFlag_) {
		vector<int> bcPts;
		//first sweep through points to find all boundary points
		for (size_t i = 0; i < points_.size(); i++) {
			if (bcFlags_[i] == 0) {
				for (int j = 0; j < adjacency[i].size(); j++) {
					if (bcFlags_[adjacency[i].at(j)] == 2) {
						for (int k = 0; k < adjacency[adjacency[i].at(j)].size(); k++) {
							auto it = std::find(adjacency[i].begin(), adjacency[i].end(), adjacency[adjacency[i].at(j)].at(k));
							if (it == adjacency[i].end()) {
								adjacency[i].push_back(adjacency[adjacency[i].at(j)].at(k));
							}
						}
					}
				}
			}
		}

		vector<int> remove_pts;
		int count = 0;
		//std::vector<int>::iterator new_end;
		for(size_t i = 0; i < points_.size(); i++){
			if(bcFlags_[i] == 0){
				for(int j = 0; j < adjacency[i].size(); j++){
					if(bcFlags_[adjacency[i].at(j)] == 2){
						remove_pts.push_back(adjacency[i].at(j));
					}
				}

				for(int k = 0; k < remove_pts.size(); k++){
					auto it = std::find(adjacency[i].begin(), adjacency[i].end(),remove_pts[k]);
					adjacency[i].erase(it);
				}

				remove_pts.clear();

			}
			else{
				for(int j = 0; j < adjacency[i].size(); j++){
					if(bcFlags_[adjacency[i].at(j)] == 0){
						remove_pts.push_back(adjacency[i].at(j));
					}

				}
				
				for(int k = 0; k < remove_pts.size(); k++){
					auto it = std::find(adjacency[i].begin(), adjacency[i].end(),remove_pts[k]);
					adjacency[i].erase(it);
				}


				
				remove_pts.clear();

				
			}
		}
	}

	reverse_cuthill_mckee_ordering(adjacency, order,bcFlags_);

	vector<Point> newPoints = points_;
	vector<Point> newNormVecs = normalVecs_;
	Eigen::VectorXd newSource = source_;
	vector<int> newBCFlag = bcFlags_;
	vector<int> oldToRCMPtrs(points_.size());

	for (size_t i = 0; i < points_.size(); i++) {
		newPoints[i] = points_[order[i]];
		newBCFlag[i] = bcFlags_[order[i]];
		newSource(i) = source_.coeff(order[i]);
		newNormVecs[i] = normalVecs_[order[i]];
		oldToRCMPtrs[order[i]] = i;
	}
	points_ = newPoints;
	source_ = newSource;
	bcFlags_ = newBCFlag;
	normalVecs_ = newNormVecs;

	vector<int> newBCPoints;
	vector<double> newBCVals;

	for (size_t i = 0; i < boundaries_.size(); i++) {
		newBCPoints.clear();
		newBCVals.clear();
		for (size_t j = 0; j < boundaries_[i].bcPoints.size(); j++) {
			newBCPoints.push_back(oldToRCMPtrs[boundaries_[i].bcPoints.at(j)]);
			newBCVals.push_back(boundaries_[i].values.at(j));
		}
		boundaries_[i].bcPoints = newBCPoints;
		boundaries_[i].values = newBCVals;
	}

}

void Grid::build_interior_mat(){
	vector<Eigen::Triplet<double>> tripletList;
	vector<Eigen::Triplet<double>> tripletList2;
	int count = 0;
	for (int i = 0; i < laplaceMatSize_; i++ ){
		if(bcFlags_[i] == 0){
			tripletList.push_back(Eigen::Triplet<double>(count,i,1));
			tripletList2.push_back(Eigen::Triplet<double>(i,count,1));
			count++;
		}
	}

	if(regularizeFlag_){
		tripletList.push_back(Eigen::Triplet<double>(count,laplaceMatSize_,1));
		tripletList2.push_back(Eigen::Triplet<double>(laplaceMatSize_,count,1));
		count++;
	}





	restrict_mat->setFromTriplets(tripletList.begin(), tripletList.end());
	restrict_mat->makeCompressed();
	prolong_mat->setFromTriplets(tripletList2.begin(), tripletList2.end());
	prolong_mat->makeCompressed();


	*(interior_mat) = *(restrict_mat) * *(laplaceMat_) * *(prolong_mat);

}

int Grid::getSize() {
	return laplaceMatSize_;
}
int Grid::getStencilSize() {
	return properties_.stencilSize;
}
int Grid::getPolyDeg() {
	return properties_.polyDeg;
}