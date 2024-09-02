#include <fstream>
#include "file_io.hpp"

//file_io.cpp - чтение матрицы из файла и вывод матрицы в файл.
//чтение матрицы с файла
Matrix readMatrixFromFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    int rows, cols;
    inputFile >> rows >> cols;
    Matrix matrix = createMatrix(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            inputFile >> matrix.data[i][j];
        }
    }

    inputFile.close();
    return matrix;
}

//запись матрицы в файл
void writeMatrixToFile(const Matrix& matrix, const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    outputFile << matrix.rows << ' ' << matrix.cols << '\n';
    for (int i = 0; i < matrix.rows; ++i) {
        for (int j = 0; j < matrix.cols; ++j) {
            outputFile << matrix.data[i][j] << ' ';
        }
        outputFile << '\n';
    }

    outputFile.close();
}