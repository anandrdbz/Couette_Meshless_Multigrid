#include "general_computation_functions.h"
using std::cout;
using std::endl;
double distance(Point refPoint, Point queryPoint) {
	return std::sqrt(std::pow(std::get<0>(refPoint) - std::get<0>(queryPoint), 2) + std::pow(std::get<1>(refPoint) - std::get<1>(queryPoint), 2));
}

std::pair<double, double> minMaxCoord(vector<Point> points_, char coord){
	double min;
	double max;
	double tempCoord;
	if (coord == 'x') {
		min = std::get<0>(points_[0]);
		max = min;
	}
	else if (coord == 'y') {
		min = std::get<1>(points_[0]);
		max = min;
	}
	for (size_t i = 0; i < points_.size(); i++) {
		if (coord == 'x') {
			tempCoord = std::get<0>(points_[i]);
		}
		else if (coord == 'y') {
			tempCoord = std::get<1>(points_[i]);
		}
		if (tempCoord > max) {
			max = tempCoord;
		}
		else if (tempCoord < min) {
			min = tempCoord;
		}
	}
	return std::pair<double, double>(min, max);
}
Point vec_from_pts(Point point1, Point point2) {
	double xdiff = std::get<0>(point2) - std::get<0>(point1);
	double ydiff = std::get<1>(point2) - std::get<1>(point1);
	double zdiff = std::get<2>(point2) - std::get<2>(point1);

	//cout << "xdiff: " << xdiff << " ydiff: " << ydiff << endl;
	return std::make_tuple(xdiff, ydiff, zdiff);
}
Point midpoint(Point point1, Point point2) {
	double xdiff = std::get<0>(point2) + std::get<0>(point1);
	double ydiff = std::get<1>(point2) + std::get<1>(point1);
	double zdiff = std::get<2>(point2) + std::get<2>(point1);

	//cout << "xdiff: " << xdiff << " ydiff: " << ydiff << endl;
	return std::make_tuple(xdiff/2, ydiff/2, zdiff/2);
}

Point centroid(Point point1, Point point2, Point point3) {
	double xdiff = std::get<0>(point3) + std::get<0>(point2) + std::get<0>(point1);
	double ydiff = std::get<0>(point3) + std::get<1>(point2) + std::get<1>(point1);

	return std::make_tuple(xdiff / 3, ydiff / 3, 0);
}

Point unit_normal_vec(Point vec, bool cw) {
	double x = std::get<0>(vec);
	double y = std::get<1>(vec);

	double norm = std::sqrt(x*x + y * y);
	x /= norm;
	y /= norm;
	return cw? std::make_tuple(-y, x, 0) : std::make_tuple(y, -x, 0);
}
Point avg_unit_norm_vec(Point norm1, Point norm2) {
	Point vec = norm1;
	vec = std::make_tuple(std::get<0>(norm1) + std::get<0>(norm2),
		std::get<1>(norm1) + std::get<1>(norm2), 0);

	double x = std::get<0>(vec);
	double y = std::get<1>(vec);

	double norm = std::sqrt(x*x + y * y);
	x /= norm;
	y /= norm;
	return std::make_tuple(-y, x, 0);
}
std::vector<Point> shifting_scaling(vector<Point> points_, Point evalPoint) {
	std::pair<double, double> minMaxX, minMaxY;
	minMaxX = minMaxCoord(points_, 'x');
	minMaxY = minMaxCoord(points_, 'y');
	double x, y, x_ss, y_ss, scale, minX, minY, maxX, maxY;
	minX = minMaxX.first;
	maxX = minMaxX.second;
	minY = minMaxY.first;
	maxY = minMaxY.second;
	scale = std::max(maxX - minX, maxY - minY);
	std::vector<Point> scaledPoints;
	for (size_t i = 0; i < points_.size(); i++) {
		x = std::get<0>(points_[i]);
		y = std::get<1>(points_[i]);
		x_ss = (x - minX) / scale;
		y_ss = (y - minY) / scale;
		scaledPoints.push_back(Point(x_ss, y_ss, 0));
	}
	scaledPoints.push_back(Point(scale, scale, scale));
	x = std::get<0>(evalPoint);
	y = std::get<1>(evalPoint);
	x_ss = (x - minX) / scale;
	y_ss = (y - minY) / scale;
	scaledPoints.push_back(Point(x_ss, y_ss, 0));
	return scaledPoints;
}
void cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order, vector <int> &bcFlags) {

	vector<bool> visited(order.size(), false);
	std::queue<int> queue;
	vector<int> new_order;
	int startPoint = 0;
	while(bcFlags[startPoint] != 0){
		startPoint++;
	}
	visited[startPoint] = true;
	queue.push(startPoint);
	int currPoint, adjPoint;
	while (!queue.empty()) {
		currPoint = queue.front();
		new_order.push_back(currPoint);
		queue.pop();
		for (size_t i = 0; i < adjacency[currPoint].size(); i++) {
			adjPoint = adjacency[currPoint].at(i);
			if (visited[adjPoint] == false) {
				visited[adjPoint] = true;
				queue.push(adjPoint);
			}
		}
	}

	for (int i = 0; i < order.size(); i++){
		if(bcFlags[i] == 2){
			new_order.push_back(i);
		}
	}
	order = new_order;
}
void reverse_cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order , vector<int> &bcFlags) {
	cuthill_mckee_ordering(adjacency, order, bcFlags);
	std::reverse(order.begin(), order.end());
}
