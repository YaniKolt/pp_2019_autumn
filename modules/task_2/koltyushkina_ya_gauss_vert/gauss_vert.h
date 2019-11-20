// Copyright 2019 Koltyushkina Yanina

#ifndef MODULES_TASK_2_KOLTYUSHKINA_YA_GAUSS_VERT_GAUSS_VERT_
#define MODULES_TASK_2_KOLTYUSHKINA_YA_GAUSS_VERT_GAUSS_VERT_

#include <vector>

std::vector <double> RandomMatrix(int _size);
void swop(int it, int max, int _size, std::vector <double>& mtr);
int maxind(int it, int _size, std::vector <double> mtr);
std::vector<double> PrGauss(std::vector<double> mtr, int _size);
void ObrGauss(std::vector<double> mtr, int _size, std::vector<double>& dres);
int Proverka(std::vector<double> mtr, int _size);
bool NullStr(std::vector<double> mtr, int str, int _size);
void All(std::vector<double> mtr, int _size, std::vector<double>& dres);
#endif  // MODULES_TASK_2_KOLTYUSHKINA_YAGAUSS_VERT_GAUSS_VERT_
