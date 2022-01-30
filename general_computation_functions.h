#ifndef GENERAL_COMPUTATION_H
#define GENERAL_COMPUTATION_H
#include <tuple>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <iostream>
typedef std::tuple<double, double, double> Point;
using std::vector;
double distance(Point refPoint, Point queryPoint);
std::pair<double, double> minMaxCoord(vector<Point> pointIDs, char coord);
std::vector<std::tuple<double, double, double>> shifting_scaling(vector<Point> points, Point evalPoint);
Point vec_from_pts(Point point1, Point point2);
Point midpoint(Point point1, Point point2);
Point centroid(Point p1, Point p2, Point p3);
Point unit_normal_vec(Point vec, bool cw);
Point avg_unit_norm_vec(Point norm1, Point norm2);
void cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order, vector<int> &bcFlags);
void reverse_cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order,vector<int> &bcFlags);

#endif