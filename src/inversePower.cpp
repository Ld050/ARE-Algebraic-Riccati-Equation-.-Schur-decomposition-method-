#include "matrix.hpp"
#include "qr.hpp"
#include "gauss.hpp"
#include <cstdlib> // для функции rand() и srand()
#include <ctime>   // для функции time()
#include <iostream>
#include <random>

//Метод обратных степенных итераций (A - матрица, у которой ищем собственные вектора. eigenvalue -  собственное значение, по которому ищем собственный вектор, maxIter-максимальное число итераций, Eigenvector - матрица для записи найденного вектора)
void inversePower(Matrix& A, float eigenvalue, int maxIter, Matrix& eigenvector) {
    int n = A.rows;
    Matrix Y, B;

    // Создаем матрицу X как некоторую матрицу * собственное значение
    Matrix X = createIdentityMatrix(n);
    X = scalar_multiply(X, eigenvalue);

    // Создаем матрицу Y как разницу между A и собственным значением * I
    Y = subtraction(A, X);

    deleteMatrix(X); // Освобождаем память

    // Инициализируем случайный начальный вектор B (столбец)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(1, 10);
    B = createMatrix(n, 1);
    for (int i = 0; i < n; i++) {
        B.data[i][0] = distr(gen); // случайные значения от 1 до 10
    }

    // Главный цикл обратных итераций
    Matrix T = createMatrix(n, 1);
    for (int i = 0; i < maxIter; i++) {
        Matrix newT = gauss(Y, B);  // решаем систему линейных уравнений

        deleteMatrix(T);  // удаляем старую матрицу T
        T = newT;

        // Нормировка
        double c = Norm(T); 
        T = scalar_multiply(T, 1 / c);
        B = T;
    }

    // Копируем значения из T в результирующий вектор eigenvector
    for (int i = 0; i < n; i++) {
        eigenvector.data[i][0] = B.data[i][0];
    }


}



//Сортировка собственных значений по возрастанию
void sortEigen(float* arr, int size){
    for (int i = 0; i < size; i++) {
        bool flag = true;
        for (int j = 0; j < size - (i + 1); j++) { 
            if (arr[j] > arr[j + 1]) {
            flag = false;
            std::swap(arr[j], arr[j + 1]);
            }
        }
        if (flag) {
         break;
        }
    }
}

//Объединение собственных векторов в матрицу. 
void UnionVecs(Matrix& H, Matrix& H_pure, Matrix& Result){
    Matrix Vec = createMatrix(H.rows, 1);
    float* eigenvalues = new float[H.rows];
    for(int k = 0; k < H.rows; ++k){
        eigenvalues[k] = H.data[k][k];
    }

    sortEigen(eigenvalues, H.rows);

    for(int i = 0; i < H.rows; ++i){
        inversePower(H_pure, eigenvalues[i], 100, Vec);
        for(int j = 0; j < Result.rows; ++j){
            Result.data[j][i] = Vec.data[j][0];
        }
    }

    delete[] eigenvalues;  // освобождаем eigenvalues
    deleteMatrix(Vec);     // освобождаем Vec
}