#ifndef GAUSS_H
#define GAUSS_H

#include "matrix.hpp"
#include <string>

Matrix gauss(Matrix A, Matrix B);
void backSubstitute(const Matrix &A, Matrix &x, const Matrix &b);
Matrix invertMatrix(Matrix A);
void invertSubstitute(Matrix& A, Matrix& I);

#endif // FILE_IO_H