#ifndef FILE_READING_FUNCTIONS_H
#define FILE_READING_FUNCTIONS_H
#include <vector>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
std::vector<std::tuple<double, double, double>> pointsFromMshFile(const char* fname);
std::vector<std::tuple<double, double, double>> pointsFromTxts(const char* fname);
std::vector<std::pair<int, int>> boundPtsConnFromMsh(const char* fname, const std::vector<int> & bcFlags);
std::vector<int> orderFromTxt(const char* fname, int nv);
void writeVectorToTxt(std::vector<double> vec, const char* filename);

#endif