#ifndef QR_HPP
#define QR_HPP

#include "matrix.hpp"
#include <utility>

Matrix createIdentityMatrix(int size);
double vectorNorm(const Matrix& v);
void qr_algorithm(Matrix& A, Matrix& U, int iterations);
void qr_householder_decomp(Matrix& A, Matrix& Q, Matrix& H);
void inversePower(Matrix& A, float eigenvalue, int maxIter, Matrix& eigenvector);
double Norm(const Matrix& A);
#endif