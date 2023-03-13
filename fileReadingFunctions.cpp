#include "fileReadingFunctions.h"
//stolen from Shantanu's code
using std::cout;
using std::endl;
std::vector<std::tuple<double, double, double>> pointsFromMshFile(const char* fname) {
	std::vector<std::tuple<double, double, double>> points;
	FILE *file;
	int itemp, nv;
	double xtemp, ytemp, ztemp;
	char temp[50];
	file = fopen(fname, "r");
	while (true)
	{
		fscanf(file, "%s ", temp);
		if (strcmp(temp, "$Nodes") == 0)
		{
			break;
		}
	}
	fscanf(file, "%i ", &nv);
	for (int iv = 0; iv < nv; iv++)
	{
		fscanf(file, "%i ", &itemp);  //vertex number
		fscanf(file, "%lf ", &xtemp); //x co-ordinate
		fscanf(file, "%lf ", &ytemp); //y co-ordinate
		fscanf(file, "%lf ", &ztemp); //z co-ordinate
		//std::cout << xtemp << " " << ytemp << " " << ztemp << std::endl;
		points.push_back(std::tuple<double, double, double>(xtemp, ytemp, ztemp));
		char* buf = new char [50];
		fgets(buf,50,file);
		delete[] buf;
	}
	return points;
}
std::vector<std::tuple<double, double, double>> pointsFromTxts(const char* fname) {
	std::vector<std::tuple<double, double, double>> points;
	FILE *file;
	int nv;
	double xtemp, ytemp, ztemp;
	char temp[50];
	file = fopen(fname, "r");
	while (true)
	{
		fscanf(file, "%s ", temp);
		if (strcmp(temp, "$Nodes") == 0)
		{
			break;
		}
	}
	fscanf(file, "%i ", &nv);
	for (int iv = 0; iv < nv; iv++)
	{
		fscanf(file, "%lf ", &xtemp); //x co-ordinate
		fscanf(file, "%lf ", &ytemp); //y co-ordinate
		fscanf(file, "%lf ", &ztemp); //z co-ordinate
		points.push_back(std::tuple<double, double, double>(xtemp, ytemp, ztemp));
	}
	return points;
}
std::vector<int> orderFromTxt(const char* fname, int nv) {
	std::vector<int> order;
	FILE *file;
	int itemp;
	char temp[50];
	file = fopen(fname, "r");
	for (int iv = 0; iv < nv; iv++)
	{
		fscanf(file, "%i ", &itemp); //x co-ordinate
	}
	return order;
}
void writeVectorToTxt(std::vector<double> vec, const char* filename)
{
	std::ofstream file;
	file.open(filename);
	//f << fixed << setprecision(2) << endl;
	for (size_t i = 0; i < vec.size(); i++) {
		file << vec[i] << "\n";
	}
	file.close();
}
std::vector<std::pair<int, int>> boundPtsConnFromMsh(const char* fname, const std::vector<int> & bcFlags) {
	FILE *file;
	int num_elems;
	int elem_type, curr_node, temp1, boundIndex, num_bcPoints_temp;
	std::vector<int> currBounds;
	char temp[50];
	std::vector<std::pair<int, int>> boundPtsConn(bcFlags.size(), std::make_pair(-1, -1));
	file = fopen(fname, "r");
	while (true)
	{
		fscanf(file, "%s ", temp);
		if (strcmp(temp, "$Elements") == 0)
		{
			break;
		}
	}
	fscanf(file, "%i ", &num_elems);
	for (int iv = 0; iv < num_elems; iv++)
	{
		fscanf(file, "%i ", &temp1); //Elem number
		fscanf(file, "%i ", &elem_type); //Elem type
		//for triangles
		if (elem_type == 2) {
			currBounds.clear();
			num_bcPoints_temp = 0;
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
			for (int j = 0; j < 3; j++) {
				fscanf(file, "%i ", &curr_node);
				curr_node--;

				if (bcFlags.at(curr_node) != 0) {
					currBounds.push_back(curr_node);
					num_bcPoints_temp++;
				}
			}
			//update the bc points connectivity if there are more than 2 bc points in the element
			if (num_bcPoints_temp >= 2) {
				for (int b = 0; b < num_bcPoints_temp; b++) {
					boundIndex = currBounds[b];
					for (int j = 0; j < currBounds.size(); j++) {
						if (currBounds[j] != boundIndex && boundPtsConn[boundIndex].first == -1) {
							boundPtsConn[boundIndex].first = currBounds[j];
						}
						else if (currBounds[j] != boundIndex && boundPtsConn[boundIndex].second == -1) {
							boundPtsConn[boundIndex].second = currBounds[j];
						}
					}
				}
			}
		}
		//for line elements
		else if (elem_type == 1) {
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
		}
		//for points
		else if (elem_type == 15) {
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
			fscanf(file, "%i ", &temp1);
		}

	}
	return boundPtsConn;
}
