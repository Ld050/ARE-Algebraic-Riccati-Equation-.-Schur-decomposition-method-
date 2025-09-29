#include <iostream>
#include <vector>
#include "matrix.hpp"

void backSubstitute(const Matrix &A, Matrix &x, const Matrix &b);
void invertSubstitute(Matrix& A, Matrix& I);

Matrix gauss(Matrix A, Matrix B) {
    int n = A.rows;
    Matrix x = createMatrix(n, 1); // Это будет решение xi

    for (int i = 0; i < n; i++) {
        int max_row = i;
        // Выделение максимального элемента
        for (int k = i + 1; k < n; k++) {
            if (std::abs(A.data[k][i]) > std::abs(A.data[max_row][i])) {
                max_row = k;
            }
        }

        // Перестанвока
        for (int k = i; k < n; k++) {
            std::swap(A.data[max_row][k], A.data[i][k]);
        }
        std::swap(B.data[max_row][0], B.data[i][0]);

        // Прямой ход
        for (int k = i + 1; k < n; k++) {
            long double factor = A.data[k][i] / A.data[i][i];
            for (int j = i; j < n; j++) {
                A.data[k][j] -= A.data[i][j] * factor;
            }
            B.data[k][0] -= B.data[i][0] * factor;
        }
    }

    backSubstitute(A, x, B);
    return x;
}
//обратный ход - поднимаемся с последнего уравнения к первому 
void backSubstitute(const Matrix &A, Matrix &x, const Matrix &b) {
    int n = A.rows;
    for (int i = n - 1; i >= 0; i--) {
        x.data[i][0] = b.data[i][0];
        for (int j = i + 1; j < n; j++) {
            x.data[i][0] -= A.data[i][j] * x.data[j][0];
        }
        x.data[i][0] /= A.data[i][i];
    }
}

Matrix invertMatrix(Matrix A) {
    int n = A.rows;
    Matrix I = createMatrix(n, n); 

    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                I.data[i][j] = 1.0;
            } else {
                I.data[i][j] = 0.0;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        int max_row = i;
        // Находим строку с максимальным элементом в текущем столбце
        for (int k = i + 1; k < n; k++) {
            if (std::abs(A.data[k][i]) > std::abs(A.data[max_row][i])) {
                max_row = k;
            }
        }

        // Перестановка строк
        for (int k = 0; k < n; k++) {
            std::swap(A.data[max_row][k], A.data[i][k]);
            std::swap(I.data[max_row][k], I.data[i][k]);
        }

        // Прямой ход
        for (int k = i + 1; k < n; k++) {
            long double factor = A.data[k][i] / A.data[i][i];
            for (int j = i; j < n; j++) {
                A.data[k][j] -= A.data[i][j] * factor;
            }
            for (int j = 0; j < n; j++) {
                I.data[k][j] -= I.data[i][j] * factor;
            }
        }
    }

    invertSubstitute(A, I);
    return I;
}

void invertSubstitute(Matrix& A, Matrix& I) {
    int n = A.rows;
    for (int i = n - 1; i >= 0; i--) {
        long double factor = A.data[i][i];
        for (int j = 0; j < n; j++) {
            I.data[i][j] /= factor;
            A.data[i][j] /= factor;
        }
        for (int k = i - 1; k >= 0; k--) {
            factor = A.data[k][i];
            for (int j = 0; j < n; j++) {
                I.data[k][j] -= I.data[i][j] * factor;
                A.data[k][j] -= A.data[i][j] * factor;
            }
        }
    }
}
