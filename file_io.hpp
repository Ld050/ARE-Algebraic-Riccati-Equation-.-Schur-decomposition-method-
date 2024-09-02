#ifndef FILE_IO_H
#define FILE_IO_H

#include "matrix.hpp"
#include <string>

Matrix readMatrixFromFile(const std::string& filename);
void writeMatrixToFile(const Matrix& matrix, const std::string& filename);

#endif // FILE_IO_H