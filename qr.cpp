#include "qr.hpp"
#include <cmath>
#include <utility>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>


//qr.cpp - прописан qr алгоритм на преобразованиях Хаусхолдера.
//Создание единичной матрицы (элементы на главной диагонали = 1)
Matrix createIdentityMatrix(int size) {
    Matrix I = createMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            I.data[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    return I;
}

// Функция для нахождения евклидовой нормы вектора
double vectorNorm(const Matrix& v) {
    double sum = 0.0;
    for (int i = 0; i < v.rows; ++i) {
        sum += v.data[i][0] * v.data[i][0];
    }
    return sqrt(sum);
}

//норма матричных значений
double Norm(const Matrix& A){
    float sum = 0;
    for(int i = 0; i < A.rows; ++i){
        for(int j = 0; j < A.cols; ++j){
            sum = sum + A.data[i][j]*A.data[i][j];
        }
    }
    sum = sqrt(sum);
    return sum;
}

//QR разложение матрицы на преобразованиях Хаусхолдера.
void qr_householder_decomp(Matrix& A, Matrix& Q, Matrix& H) {
    Q = createIdentityMatrix(A.rows);
    double col_quad_sum = 0.0;
    Matrix v = createMatrix(A.rows, 1);
    Matrix v_product, v_scalar;
    //запускаем основной цикл для H_n преобразований. 
    for (int i_main = 0; i_main < (A.cols - 1); ++i_main) {
        for (int i = 0; i < A.rows; ++i) {
            v.data[i][0] = 0.0;
        }

        col_quad_sum = 0.0;
        for (int i = i_main; i < A.rows; ++i) {
            col_quad_sum += A.data[i][i_main] * A.data[i][i_main];
        }
        v.data[i_main][0] = A.data[i_main][i_main] + copysign(sqrt(col_quad_sum), A.data[i_main][i_main]);

        for (int i = (i_main + 1); i < A.rows; ++i) {
            v.data[i][0] = A.data[i][i_main];
        }

        // Создаем временные матрицы для v и освобождаем ранее созданные сразу после использования
        Matrix temp = transposeMatrix(v);
        v_product = multiplyMatrix(v, temp);
        v_scalar = multiplyMatrix(temp, v);
        deleteMatrix(temp); // Освобождаем память, использованную для temp (по сути транспонированной v).
        //вычисляем v
        double scalar = 2.0 / v_scalar.data[0][0];
        for (int i = 0; i < v_product.rows; ++i) {
            for (int j = 0; j < v_product.cols; ++j) {
                v_product.data[i][j] *= scalar;
            }
        }
        //вычисляем H
        H = createIdentityMatrix(A.rows);
        for (int i = 0; i < H.rows; ++i) {
            for (int j = 0; j < H.cols; ++j) {
                H.data[i][j] -= v_product.data[i][j];
            }
        }

        // Создаем временные матрицы для умножения Q и H и освобождаем ранее созданные сразу после использования
        temp = multiplyMatrix(Q, H);
        deleteMatrix(Q);  // Удаляем старую матрицу Q
        Q = temp;         // Присваиваем новую матрицу Q

        temp = multiplyMatrix(H, A);
        deleteMatrix(A);  // Удаляем старую матрицу A
        A = temp;         // Присваиваем новую матрицу A

        deleteMatrix(H);  // Удаляем H в конце каждой итерации
        deleteMatrix(v_product);
        deleteMatrix(v_scalar);
    }

    deleteMatrix(v);
}
void qr_algorithm(Matrix& A, Matrix& U, int iterations) {
    Matrix E = createMatrix(A.rows, A.cols); // Инициализируем матрицу E
    U = createIdentityMatrix(A.rows);
    Matrix Q, H;

    for (int i = 0; i < iterations; ++i) {
        E = createIdentityMatrix(A.rows); 
        for (int j = 0; j < E.rows; ++j) {
            E.data[j][j] *= A.data[A.rows - 1][A.cols - 1]; 
        }

        Matrix temp = subtraction(A, E); 

        qr_householder_decomp(temp, Q, H); // QR разложение для temp

        U = multiplyMatrix(U, Q); //Должна быть матрицой собственных векторов, которая накапливает в себе произведение Q. Но работает не для всех матриц.
        //Стабильный метод для поиска собственных векторов для матриц, у которых размерность меньше чем 5x5 - обратные степенные итерации.

        Matrix temp2 = multiplyMatrix(temp, Q); // A_new = temp * Q_new
        A = addMatrix(temp2, E); // A = A_new + E
        deleteMatrix(temp);
        deleteMatrix(temp2);
        deleteMatrix(Q); // освобождаем Q после завершения итерации
        deleteMatrix(H); // освобождаем H после завершения итерации
    }
    deleteMatrix(E); // Почистить выделенную память
}
































































































/*
// QR-разложение методом Хаусхолдера
// QR-разложение методом Хаусхолдера
std::pair<Matrix, Matrix> qrDecomposition(const Matrix& A) {
    int m = A.rows;
    int n = A.cols;
    Matrix Q = createIdentityMatrix(m);
    Matrix R = createMatrix(m, n);

    // Копирование исходной матрицы A в R
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            R.data[i][j] = A.data[i][j];
        }
    }

    for (int k = 0; k < n && k < m - 1; ++k) {
        // Создаем вектор x
        Matrix x = createMatrix(m - k, 1);
        for (int i = k; i < m; ++i) {
            x.data[i - k][0] = R.data[i][k];
        }

        // Норма вектора x
        double norm_x = vectorNorm(x);

        // Создаем вектор u
        Matrix u = createMatrix(m - k, 1);
        for (int i = 0; i < m - k; ++i) {
            u.data[i][0] = x.data[i][0];
        }
        u.data[0][0] += (x.data[0][0] >= 0) ? norm_x : -norm_x;

        // Норма вектора u
        double norm_u = vectorNorm(u);

        if (norm_u != 0.0) {
            for (int i = 0; i < m - k; ++i) {
                u.data[i][0] /= norm_u;
            }
        }

        // Создаем матрицу Hk
        Matrix Hk = createIdentityMatrix(m - k);
        Matrix uuT = createMatrix(m - k, m - k);
        for (int i = 0; i < m - k; ++i) {
            for (int j = 0; j < m - k; ++j) {
                uuT.data[i][j] = u.data[i][0] * u.data[j][0];
            }
        }

        for (int i = 0; i < m - k; ++i) {
            for (int j = 0; j < m - k; ++j) {
                Hk.data[i][j] -= 2.0 * uuT.data[i][j];
            }
        }

        // Расширяем Hk до размера m x m
        Matrix H = createIdentityMatrix(m);
        for (int i = k; i < m; ++i) {
            for (int j = k; j < m; ++j) {
                H.data[i][j] = Hk.data[i - k][j - k];
            }
        }

        // Обновляем R и Q
        Matrix R_temp = multiplyMatrix(H, R);
        deleteMatrix(R);
        R = R_temp;

        Matrix Q_temp = multiplyMatrix(Q, H);
        deleteMatrix(Q);
        Q = Q_temp;

        deleteMatrix(x);
        deleteMatrix(u);
        deleteMatrix(Hk);
        deleteMatrix(H);
        deleteMatrix(uuT);
    }

    return std::make_pair(Q, R);
}

std::pair<Matrix, Matrix> qrAlgorithm(const Matrix& A) {
    Matrix H = A;
    Matrix U = createIdentityMatrix(A.rows); // Начальная матрица U (единичная)

    for (int i = 0; i < 10000; ++i) { // Итерируем, чтобы приблизить к форме Шура
        auto [Q, R] = qrDecomposition(H);
        
        Matrix H_temp = multiplyMatrix(R, Q);
        deleteMatrix(H);
        H = H_temp;

        Matrix U_temp = multiplyMatrix(U, Q);
        deleteMatrix(U);
        U = U_temp;

        deleteMatrix(Q);
        deleteMatrix(R);
    }
    return std::make_pair(U, H);
}
*/