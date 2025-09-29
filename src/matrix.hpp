#ifndef MATRIX_H
#define MATRIX_H

struct Matrix {
    int rows;
    int cols;
    double** data;
};


Matrix createMatrix(int rows, int cols);
void deleteMatrix(Matrix& matrix);
Matrix addMatrix(const Matrix& A, const Matrix& B);
Matrix subtraction(const Matrix& A, const Matrix& B);
Matrix multiplyMatrix(const Matrix& A, const Matrix& B);
Matrix transposeMatrix(const Matrix& A);
Matrix destroy_nills(const Matrix& A, float precision);
Matrix inverse(const Matrix& A);
Matrix scalar_multiply(const Matrix& A, double scalar);
void printMatrix(const Matrix& matrix);
void printDiag(const Matrix& A);
double get_determinant(const Matrix& A);
void constructBlocks(const Matrix& A, Matrix& B_1, Matrix& B_2);

#endif // MATRIX_H