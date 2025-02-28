#ifndef TASK_H
#define TASK_H
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>

int gaussJordanInverse(std::vector<std::vector<double>>& A,
    std::vector<std::vector<double>>& inv,std::vector<int> &columnOrder, std::vector<int> &undo);
#endif
