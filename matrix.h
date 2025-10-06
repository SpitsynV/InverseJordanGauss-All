#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>



int readMatrixFromFile(const std::string& filename, std::vector<std::vector<double>>& A, const int n);
void printMatrix(const std::vector<std::vector<double>>& A, const int m);
void printVector(const std::vector<double>& vec, int m);
inline double f(const int k, const int n,const  int i,const  int j);
void initializeMatrix(std::vector<std::vector<double>>& A, const int k, const int n);

#endif
