#include "testing_functions.hpp"

double calc_l1_error(Grid* grid, bool neumannFlag, int k1, int k2) {
	vector<Point> points = grid->points_;
	Eigen::VectorXd actual(points.size());
	double x, y;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		actual(i) = neumannFlag ? std::cos(k1 * pi * x) * std::cos(k2 * pi * y) : std::sin(k1 * pi * x) * std::sin(k2 * pi * y);
	}
	double solMean = 0;
	double manufacturedMean = 0;
	if (!neumannFlag) {
		return (*grid->values_ - actual).lpNorm<1>() / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		solMean += grid->values_->coeff(i) / grid->laplaceMatSize_;
		manufacturedMean += actual.coeff(i) / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		grid->values_->coeffRef(i) += (manufacturedMean - solMean);
	}

	double error = 0;
	for (int i = 0; i < points.size(); i++) {
		error += std::abs(grid->values_->coeff(i) - actual.coeff(i));
	}
	error /= points.size();
	return error;
}
double calc_l1_error_circle(Grid* grid, bool neumannFlag, int k) {
	vector<Point> points = grid->points_;
	Eigen::VectorXd actual(points.size());
	double x, y, rstar;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		x -= 0.5;
		y -= 0.5;
		rstar = (std::sqrt(x * x + y * y) - 0.25) / (0.5 - 0.25);
		actual(i) = std::sin(rstar * pi * k);
	}
	double solMean = 0;
	double manufacturedMean = 0;
	if (!neumannFlag) {
		return (*grid->values_ - actual).lpNorm<1>() / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		solMean += grid->values_->coeff(i) / grid->laplaceMatSize_;
		manufacturedMean += actual.coeff(i) / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		grid->values_->coeffRef(i) += (manufacturedMean - solMean);
	}

	double error = 0;
	for (int i = 0; i < points.size(); i++) {
		error += std::abs(grid->values_->coeff(i) - actual.coeff(i));
	}
	error /= points.size();
	return error;
}
Grid* genGmshGridDirichlet(string geomtype, const char* filename, GridProperties props, std::string filetype, int k1, int k2) {
	vector<std::tuple<double, double, double>> points;
	if (filetype.compare("txt") == 0) {
		points = pointsFromTxts(filename);
	}
	else {
		points = pointsFromMshFile(filename);
	}
	double x, y;
	vector<int> bPts, bPts_inner;
	vector<double> bValues, bValues_inner;

	Eigen::VectorXd source(points.size());
	if (geomtype.compare("square") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			source(i) = -(k1 * k1 + k2 * k2) * pi * pi * std::sin(k1 * pi * x) * std::sin(k2 * pi * y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0.0);
			}
		}
	}
	else if (geomtype.compare("square_with_circle") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			source(i) = -(k1 * k1 + k2 * k2) * pi * pi * std::sin(k1 * pi * x) * std::sin(k1 * pi * y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0);
			}
			else if (std::abs(0.0625 - (x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5)) <= std::pow(10, -10)) {
				bPts_inner.push_back(i);
				bValues_inner.push_back(std::sin(k1 * pi * x) * std::sin(k1 * pi * y));
			}
		}
	}
	else if (geomtype.compare("concentric_circles") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			x -= 0.5;
			y -= 0.5;
			double sum = 0;
			double r = std::sqrt(x * x + y * y);
			double rstar = (r - 0.25) / (0.5 - 0.25);
			//cout << I << endl;
			sum += -pi * k1 * k1 * pi * std::sin(pi * k1 * rstar) * std::pow(4 * x * std::pow(x * x + y * y, -0.5), 2)
				+ pi * k1 * std::cos(pi * k1 * rstar) * 4 * (std::pow(x * x + y * y, -0.5) + 2 * x * x * -0.5 * std::pow(x * x + y * y, -1.5));
			sum += -pi * k1 * k1 * pi * std::sin(pi * k1 * rstar) * std::pow(4 * y * std::pow(x * x + y * y, -0.5), 2)
				+ pi * k1 * std::cos(pi * k1 * rstar) * 4 * (std::pow(x * x + y * y, -0.5) + 2 * y * y * -0.5 * std::pow(x * x + y * y, -1.5));
			x += 0.5;
			y += 0.5;
			source(i) = sum;
			if (std::abs(0.25 - (x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5)) <= std::pow(10, -10)) {
				bPts.push_back(i);
				bValues.push_back(0.0);

			}
			else if (std::abs(0.0625 - (x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5)) <= std::pow(10, -10)) {
				bPts_inner.push_back(i);
				bValues_inner.push_back(0.0);

			}
		}
	}
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 1;
	boundary.values = bValues;
	vector<Boundary> bcs;
	bcs.push_back(boundary);
	if (geomtype.compare("square") != 0) {
		Boundary inner;
		inner.bcPoints = bPts_inner;
		inner.type = 1;
		inner.values = bValues_inner;
		bcs.push_back(inner);
	}

	Grid* grid = new Grid(points, bcs, props, source);
	grid->implicitFlag_ = false;
	grid->setBCFlag(0, std::string("dirichlet"), bValues);
	if (geomtype.compare("square") != 0) {
		grid->setBCFlag(1, std::string("dirichlet"), bValues_inner);
	}
	grid->rcm_order_points();
	grid->build_laplacian();
	return grid;
}

Grid* genGmshGridNeumann(string geomtype, const char* filename, GridProperties props, std::string filetype, int k1, int k2, std::string coarse) {
	vector<std::tuple<double, double, double>> points;
	if (filetype.compare("txt") == 0) {
		points = pointsFromTxts(filename);
	}
	else {
		points = pointsFromMshFile(filename);
	}
	double x, y;
	vector<int> bPts, bPts_inner;
	vector<double> bValues, bValues_inner;
	Eigen::VectorXd actual(points.size() + 1);
	Eigen::VectorXd source(points.size() + 1);
	if (geomtype.compare("square") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			actual(i) = std::cos(k1 * pi * x) * std::cos(k2 * pi * y);
			source(i) = -(k1 * k1 + k2 * k2) * pi * pi * std::cos(k1 * pi * x) * std::cos(k2 * pi * y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0.0);
			}
		}
	}
	else if (geomtype.compare("square_with_circle") == 0) {
		double normX, normY, norm;
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			source(i) = -(k1 * k1 + k2 * k2) * pi * pi * std::cos(k1 * pi * x) * std::cos(k2 * pi * y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0);
			}
			else if (std::abs(0.0625 - (x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5)) <= std::pow(10, -10)) {
				normX = x - 0.5;
				normY = y - 0.5;
				norm = std::sqrt(normX * normX + normY * normY);
				normX /= norm;
				normY /= norm;
				bPts_inner.push_back(i);
				bValues_inner.push_back(-normX * pi * k1 * std::sin(k1 * pi * x) * std::cos(k2 * pi * y)
					- normY * pi * k2 * std::cos(k1 * pi * x) * std::sin(k2 * pi * y));

			}
		}
	}
	else if (geomtype.compare("concentric_circles") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			x -= 0.5;
			y -= 0.5;
			double sum = 0;
			double r = std::sqrt(x * x + y * y);
			double rstar = (r - 0.25) / (0.5 - 0.25);
			//cout << I << endl;
			sum += -pi * k1 * k1 * pi * std::sin(pi * k1 * rstar) * std::pow(4 * x * std::pow(x * x + y * y, -0.5), 2)
				+ pi * k1 * std::cos(pi * k1 * rstar) * 4 * (std::pow(x * x + y * y, -0.5) + 2 * x * x * -0.5 * std::pow(x * x + y * y, -1.5));
			sum += -pi * k1 * k1 * pi * std::sin(pi * k1 * rstar) * std::pow(4 * y * std::pow(x * x + y * y, -0.5), 2)
				+ pi * k1 * std::cos(pi * k1 * rstar) * 4 * (std::pow(x * x + y * y, -0.5) + 2 * y * y * -0.5 * std::pow(x * x + y * y, -1.5));
			x += 0.5;
			y += 0.5;
			source(i) = sum;
			double normX, normY, norm;
			if (std::abs(0.25 - (x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5)) <= std::pow(10, -10)) {
				normX = x - 0.5;
				normY = y - 0.5;
				norm = std::sqrt(normX * normX + normY * normY);
				normX /= norm;
				normY /= norm;
				bPts.push_back(i);
				bValues.push_back(-normX * k1 * pi * std::cos(k1 * pi * rstar) / r * 4 * (x - 0.5)
					- normY * k1 * pi * std::cos(k1 * pi * rstar) / r * 4 * (y - 0.5));

			}
			else if (std::abs(0.0625 - (x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5)) <= std::pow(10, -10)) {
				normX = x - 0.5;
				normY = y - 0.5;
				norm = std::sqrt(normX * normX + normY * normY);
				normX /= norm;
				normY /= norm;

				bPts_inner.push_back(i);
				bValues_inner.push_back(normX * k1 * pi * std::cos(k1 * pi * rstar) / r * 4 * (x - 0.5)
					+ normY * k1 * pi * std::cos(k1 * pi * rstar) / r * 4 * (y - 0.5));

			}
		}
	}
	actual(actual.rows() - 1) = 0;
	source(source.rows() - 1) = 0;
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 2;
	boundary.values = bValues;
	vector<Boundary> bcs;
	bcs.push_back(boundary);
	if (geomtype.compare("square") != 0) {
		Boundary inner;
		inner.bcPoints = bPts_inner;
		inner.type = 2;
		inner.values = bValues_inner;
		bcs.push_back(inner);
	}
	Grid* grid = new Grid(points, bcs, props, source);
	grid->implicitFlag_ = true;
	grid->setBCFlag(0, std::string("neumann"), bValues);
	if (geomtype.compare("square") != 0) {
		grid->setBCFlag(1, std::string("neumann"), bValues_inner);

	}
	//grid->print_bc_values();
	grid->build_normal_vecs(filename, geomtype);
	grid->rcm_order_points();
	grid->build_deriv_normal_bound();
	grid->build_laplacian();
	grid->modify_coeff_neumann(coarse);

	grid->push_inhomog_to_rhs();

	return grid;
}
void write_temp_contour(Grid* testGrid, string directory, string extension) {
	//Write necessary things for temperature contour
	vector<Point> points = testGrid->points_;
	double x, y;
	vector<double> xv, yv;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		xv.push_back(x);
		yv.push_back(y);
	}
	string x_str = "x_";
	string y_str = "y_";
	string temp_str = "temp_";
	string cond_str = "cond_rbf_";
	string txt = ".txt";

	writeVectorToTxt(xv, (directory + x_str + extension + txt).c_str());
	writeVectorToTxt(yv, (directory + y_str + extension + txt).c_str());

	vector<double> temp;
	for (int i = 0; i < testGrid->values_->rows() - 1; i++) {
		temp.push_back(testGrid->values_->coeff(i));
	}
	writeVectorToTxt(temp, (directory + temp_str + extension + txt).c_str());
}
void write_mg_resid(Multigrid& mg, string directory, string extension) {
	string res_str = "resid_";
	string txt = ".txt";
	writeVectorToTxt(mg.residuals_, (directory + res_str + extension + txt).c_str());
}
void write_l1error_cond(Multigrid& mg, MultigridParameters params, string directory, string extension, double clock) {
	string res_str = "cond_error_";
	string txt = ".txt";
	vector<double> vec;
	Grid* grid = mg.grids_.at(mg.grids_.size() - 1).second;
	if (params.geomtype.compare("concentric_circles") == 0) {
		vec.push_back(calc_l1_error_circle(grid, grid->neumannFlag_, params.k1));
	}
	vec.push_back(calc_l1_error(grid, grid->neumannFlag_, params.k1, params.k2));
	vec.push_back(clock);
	writeVectorToTxt(vec, (directory + res_str + extension + txt).c_str());
}
void run_mg_sim(MultigridParameters params) {
	Multigrid mg;
	for (int i = 0; i < (int)(params.filenames.size()); i++) {
		if (params.neumann == true) {
			std::string coarse = (i == params.filenames.size() - 1) ? "fine" : "coarse";
			mg.addGrid(genGmshGridNeumann(params.geomtype, (params.directory + params.filenames[i]).c_str(), params.props[i], params.filetypes[i], params.k1, params.k2, coarse));
		}
		else {
			mg.addGrid(genGmshGridDirichlet(params.geomtype, (params.directory + params.filenames[i]).c_str(), params.props[i], params.filetypes[i], params.k1, params.k2));
		}
	}
	mg.buildMatrices();
	std::clock_t start = std::clock();
	for (int i = 0; i < params.num_v_cycle; i++) {
		mg.vCycle();
	}
	double vCycTime = (std::clock() - start) / (double)(CLOCKS_PER_SEC);
	write_mg_resid(mg, params.directory, params.extension);
	write_temp_contour(mg.grids_.back().second, params.directory, params.extension);

	//write cond number and l1 error
	write_l1error_cond(mg, params, params.directory, params.extension, vCycTime);
}
MultigridParameters gen_mg_param(string geom, int numGrids, int k, int poly_deg, int vcyc, bool neumann) {
	string dir;
	vector<string> msh_files;
	if (geom.compare("square") == 0) {
		msh_files = { /*"square_98.msh",*/ "square_170.msh", "square_600.msh", "square_2.5k.msh", "square_10k.msh" };
		dir = "square_test_geometries/";
	}
	else if (geom.compare("square_with_circle") == 0) {
		msh_files = {/*"square_hole_89.msh", */"square_hole_176.msh", "square_hole_640.msh", "square_hole_2532.msh", "square_hole_10197.msh", "square_hole_37943.msh", "square_hole_150214.msh" };
		dir = "square_hole_geoms/";
	}
	else if (geom.compare("concentric_circles") == 0) {
		msh_files = {/* "concentric_circles_90.msh",*/ "concentric_circles_188.msh" , "concentric_circles_650.msh" , "concentric_circles_2581.msh" , "concentric_circles_10207.msh" };
		dir = "concentric_circle_geoms/";
	}
	vector<string> filenames, filetypes;
	for (int i = 0; i < numGrids; i++) {
		filenames.push_back(msh_files.at(i));
		filetypes.push_back("msh");
		//filetypes.push_back("txt");
	}
	vector<GridProperties> props(numGrids);
	for (int i = 0; i < numGrids; i++) {
		props[i].iters = 5;
		props[i].polyDeg = (i == numGrids - 1) ? poly_deg : 3;
		props[i].omega = 1.4;
		props[i].rbfExp = 3;
		props[i].stencilSize = (i == 0) ? (int)(2.5 * (props[i].polyDeg + 1) * (props[i].polyDeg + 2) / 2)
			: (int)(2.5 * (props[i].polyDeg + 1) * (props[i].polyDeg + 2) / 2);
	}
	string bctype = (neumann) ? "neumann" : "dirichlet";
	MultigridParameters params;
	params.geomtype = geom;
	params.directory = dir;
	params.extension = std::to_string(numGrids) + "grid_" + bctype + "_L=" +
		std::to_string(poly_deg) + "_K=" + std::to_string(k);
	params.filenames = filenames;
	params.filetypes = filetypes;
	params.k1 = k;
	params.k2 = k;
	params.neumann = neumann;
	params.num_v_cycle = vcyc;
	params.props = props;
	return params;
}
void run_tests() {

	vector<MultigridParameters> params;
	
	for (int grids = 2; grids <= 4; grids++) {
		for (int k = 1; k <= 4; k++) {
			for (int l = 4; l <= 6; l++) {
				params.push_back(gen_mg_param("square", grids, k, l, 750, true));
				params.push_back(gen_mg_param("square_with_circle", grids, k, l, 750, true));
				params.push_back(gen_mg_param("concentric_circles", grids, k, l, 750, true));
				//params.push_back(gen_mg_param("square", grids, k, l, 750, true));
				//params.push_back(gen_mg_param("square_with_circle", grids, k, l, 750, true));
				//params.push_back(gen_mg_param("concentric_circles", grids, k, l, 750, true));

			}
		}
	}
	
	//params.push_back(gen_mg_param("square", 4, 1, 6, 150, true));
	for (int i = 0; i < params.size(); i++) {
		cout << params[i].geomtype << " " << params[i].k1 << " " << params[i].neumann << " " << params[i].filenames.size() << endl;
		run_mg_sim(params[i]);
	}

}

void testGmshSingleGrid() {
	GridProperties props;

	props.iters = 5;
	props.polyDeg = 6;
	props.omega = 1.4;
	props.rbfExp = 3;
	props.stencilSize = (int)(2.0 * (props.polyDeg + 1) * (props.polyDeg + 2) / 2);

	Grid* testGrid = genGmshGridDirichlet("square_with_circle", "square_hole_geoms/square_hole_10197.msh", props, "msh", 1, 1);
	//Grid* testGrid = genGmshGridNeumann("square", "square_test_geometries/square_10k.msh", props, "msh", 1, 1);
	//Grid* testGrid = genGmshGridNeumann("concentric_circles", "concentric_circle_geoms/concentric_circles_10207.msh", props, "msh", 1, 1, "fine");

	vector<double> res;
	testGrid->boundaryOp("fine");
	for (int i = 0; i < 1000; i++) {
		cout << "residual: " << testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>() << endl;
		res.push_back(testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>());
		testGrid->sor(testGrid->laplaceMat_, testGrid->values_, testGrid->source_);
	}
	testGrid->bound_eval_neumann();
	writeVectorToTxt(res, "residual.txt");
	cout << "L1 error: " << calc_l1_error(testGrid, testGrid->neumannFlag_, 1, 1) << endl;
	//cout << "L1 error: " << calc_l1_error_circle(testGrid, testGrid->neumannFlag_, 1) << endl;

	//band matrix plotting.
	/*
	Eigen::SparseMatrix<double, 1> * matrix = testGrid->laplaceMat_;
	const double* laplaceValues = matrix->valuePtr();
	//column indices of values
	const int* innerValues = matrix->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = matrix->outerIndexPtr();
	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	int numNonZeros = matrix->nonZeros();
	int inner;
	valueIdx = 0;
	innerIdx = 0;
	outerIdx = 0;

	vector<double> iv, jv;
	for (int i = 0; i < matrix->rows(); i++) {
		//diagCoeff = laplaceMat_->coeff(i, i);
		rowStartIdx = outerValues[outerIdx];
		rowEndIdx = outerValues[outerIdx + 1];
		//sum coeffs*x_i_old on the row i
		for (int j = rowStartIdx; j < rowEndIdx; j++) {
			iv.push_back(i);
			jv.push_back(innerValues[j]);
		}
		outerIdx++;
	}
	writeVectorToTxt(iv, "i.txt");
	writeVectorToTxt(jv, "j.txt");
	*/
	cout << testGrid->values_->maxCoeff() << endl;
	cout << testGrid->values_->minCoeff() << endl;
}