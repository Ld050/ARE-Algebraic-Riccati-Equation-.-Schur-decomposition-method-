#include <iostream>
#include "matrix.hpp"

//matrix.cpp - прописаны основные операции с матрицами.


//Создание матрицы
Matrix createMatrix(int rows, int cols) {
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        matrix.data[i] = new double[cols];
    }
    return matrix;
}
//Удаление матрицы
void deleteMatrix(Matrix& matrix) {
    for (int i = 0; i < matrix.rows; ++i) {
        delete[] matrix.data[i];
    }
    delete[] matrix.data;
    matrix.data = nullptr;
    matrix.rows = 0;
    matrix.cols = 0;
}
//Суммирование матриц
Matrix addMatrix(const Matrix& A, const Matrix& B) {
    Matrix result = createMatrix(A.rows, A.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }
    return result;
}
//Разность матриц
Matrix subtraction(const Matrix& A, const Matrix& B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        std::cout << "Matrix dimensions must match for subtraction." << std::endl;
        exit(1);
    }
    Matrix result = createMatrix(A.rows, B.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < B.cols; ++j) {
            result.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }
    return result;
}
//Умножение матриц
Matrix multiplyMatrix(const Matrix& A, const Matrix& B) {
    if(A.cols != B.rows){
        std::cout << "The number of columns in the first matrix should be equal to the number of rows in the second." << std::endl;
        exit(1);
    }
    Matrix result = createMatrix(A.rows, B.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < B.cols; ++j) {
            result.data[i][j] = 0;
            for (int k = 0; k < A.cols; ++k) {
                result.data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }
    return result;
}
//Умножение матрицы A на скаляр
Matrix scalar_multiply(const Matrix& A, double scalar) {
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            A.data[i][j] *= scalar;
        }
    }
    return A;
}

//Транспонирование матрицы
Matrix transposeMatrix(const Matrix& A) {
    Matrix result = createMatrix(A.cols, A.rows);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            result.data[j][i] = A.data[i][j];
        }
    }
    return result;
}
//удаление нулевых элеметов матрицы A с точностью до precision
Matrix destroy_nills(const Matrix& A, float precision) {
    Matrix B = A; 
    for (int i = 0; i < B.rows; ++i) {
        for (int j = 0; j < B.cols; ++j) {
            if (std::abs(B.data[i][j]) < precision) {
                B.data[i][j] = 0;
            }
        }
    }
    return B;
}
// поиск определителя матрицы
double get_determinant(const Matrix& A){
    double determinant = 0;

    if (A.rows != A.cols) {
        std::cout << "Matrix must be square to find determinant."<< std::endl;;
    }

    if (A.rows == 1) {
        return A.data[0][0];
    }

    if (A.rows == 2) {
        return A.data[0][0] * A.data[1][1] - A.data[0][1] * A.data[1][0];
    }

    for (int j_out = 0; j_out != A.cols; ++j_out) {
        Matrix matrix_inside = createMatrix(A.rows - 1, A.cols - 1);
        int counter = 0;

        for (int i = 1; i != A.rows; ++i) {
            for (int j = 0; j != A.cols; ++j) {
                if (j == j_out) continue;
                matrix_inside.data[i - 1][counter % (A.cols - 1)] = A.data[i][j];
                ++counter;
            }
        }

        if (j_out % 2 == 0) {
            determinant += A.data[0][j_out] * get_determinant(matrix_inside);
        }
        else {
            determinant -= A.data[0][j_out] * get_determinant(matrix_inside);
        }
        deleteMatrix(matrix_inside);
    }
    
    return determinant;
}

//Вывод диагональных элементов
void printDiag(const Matrix& A){
    for(int i = 0; i < A.rows; ++i){
        std::cout << A.data[i][i] << std::endl;
    }
}

//Поиск обратной матрицы
Matrix inverse(const Matrix& A){
    if (A.rows != A.cols) {
        std::cout << "Matrix must be square to find inverse." << std::endl;
        exit(1);
    }
    if (A.rows == 1 && A.cols == 1) {
        Matrix inverse_matrix = createMatrix(1, 1);
        inverse_matrix.data[0][0] = 1 / A.data[0][0]; 
        return inverse_matrix;
    }    
    double det = get_determinant(A);
    if (det == 0) {
        std::cout << "Matrix is singular, cannot find inverse." << std::endl;
        exit(1);
    }

    Matrix adjoint = createMatrix(A.rows, A.cols);

    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            // Вычисляем алгебраическое дополнение
            Matrix submatrix = createMatrix(A.rows - 1, A.cols - 1);
            int sub_i = 0, sub_j = 0;
            for (int row = 0; row < A.rows; ++row) {
                if (row == i) continue;
                for (int col = 0; col < A.cols; ++col) {
                    if (col == j) continue;
                    submatrix.data[sub_i][sub_j++] = A.data[row][col];
                    if (sub_j == A.cols - 1) {
                        sub_i++;
                        sub_j = 0;
                    }
                }
            }
            
            double sign = ((i + j) % 2 == 0) ? 1 : -1;
            adjoint.data[j][i] = sign * get_determinant(submatrix);
            deleteMatrix(submatrix);
        }
    }
    

    Matrix inverse_matrix = createMatrix(A.rows, A.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            inverse_matrix.data[i][j] = adjoint.data[i][j] / det;
        }
    }
    deleteMatrix(adjoint);
    return inverse_matrix;
}



//Вывод матрицы
void printMatrix(const Matrix& matrix) {
    for (int i = 0; i < matrix.rows; ++i) {
        for (int j = 0; j < matrix.cols; ++j) {
            std::cout << matrix.data[i][j] << ' ';
            }
        std::cout << std::endl;
    }
}

void constructBlocks(const Matrix& A, Matrix& B_1, Matrix& B_2){
    int n = A.rows;
    for (int i = 0; i < n/2; ++i) {
        for (int j = 0; j < n/2; ++j) {
            B_2.data[i][j] = A.data[n/2+i][j];
        }
    }      
    for (int i = 0; i < n/2; ++i) {
        for (int j = 0; j < n/2; ++j) {
            B_1.data[i][j] = A.data[i][j];
        }
    }
}