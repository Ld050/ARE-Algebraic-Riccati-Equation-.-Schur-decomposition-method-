#ifndef INVERSE_POWER_H
#define INVERSE_POWER_H

#include "matrix.hpp"
#include "qr.hpp"
#include <string>

void inversePower(Matrix& A, float eigenvalue, int maxIter, Matrix& eigenvector);
void sortEigen(float* arr, int size);
void UnionVecs(Matrix& H, Matrix& H_pure, Matrix& Result);
#endif 