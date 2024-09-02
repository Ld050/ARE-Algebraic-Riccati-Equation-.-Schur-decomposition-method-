#include <iostream>
#include <utility>
#include <algorithm>
#include "matrix.hpp"
#include "file_io.hpp"
#include "hamiltonian.hpp"
#include "qr.hpp"
#include "inversePower.hpp"
#include "gauss.hpp"
//A: n x n
//B: n x m
//Q: n x n
//R: m x m
int main() {
    //читаем матрицы
    Matrix A = readMatrixFromFile("A.txt");

    Matrix B = readMatrixFromFile("B.txt");

    Matrix Q = readMatrixFromFile("Q.txt");

    Matrix R = readMatrixFromFile("R.txt");

    Matrix H = constructHamiltonian(A, B, Q, R);
    //Matrix H = readMatrixFromFile("test.txt");

    Matrix U;
    //Создаем заготовки для блочных матриц
    Matrix B_1 = createMatrix(H.rows/2, H.cols/2);
    Matrix B_2 = createMatrix(H.rows/2, H.cols/2);


    //qr алгоритм
    Matrix eigenVec = createMatrix(H.rows, 1);
    Matrix H_eigen = H;
    Matrix Vecs = createMatrix(H.rows, H.cols);
    int iterations = 1000;
    qr_algorithm(H, U, iterations);

    //вывод собственных значений
    float* eigenvalues = new float[H.rows];
    std::cout<<"============ Eigen values ============"<<std::endl;
    for(int i = 0; i < H.rows; ++i){
        eigenvalues[i] = H.data[i][i];
    }
    sortEigen(eigenvalues, H.rows);
    for(int i = 0; i< H.rows; ++i){
        std::cout << eigenvalues[i] << std::endl;
    }

    //поиск и вывод собственных векторов
    std::cout<<"============ Eigen vectors ============"<<std::endl;
    UnionVecs(H, H_eigen, Vecs);
    printMatrix(Vecs);


    //вычисляем решение по формуле Шура.
    constructBlocks(Vecs, B_1, B_2);
    Matrix B_1_inv = createMatrix(B_1.rows, B_1.cols);
    for(int i = 0; i < B_1.rows; i++){
        for(int j = 0; j < B_1.cols; j++){
            B_1_inv.data[i][j] = B_1.data[i][j];
        }
    }
    B_1_inv = invertMatrix(B_1_inv);
    Matrix Solution = multiplyMatrix(B_2, B_1_inv);
    std::cout << "============ Solution ============"<<std::endl;
    printMatrix(Solution);


    std::cout << "============ A^TP + PA - PBR^-1B^TP + Q ============" << std::endl;
    Matrix ATP = multiplyMatrix(transposeMatrix(A),Solution);

    Matrix PA = multiplyMatrix(Solution, A);
    Matrix R_inv = createMatrix(R.rows, R.cols);
    for(int i = 0; i < R.rows; i++){
        for(int j = 0; j < R.cols; j++){
            R_inv.data[i][j] = R.data[i][j];
        }
    }
    R_inv = invertMatrix(R_inv);
    Matrix PBR_inv = multiplyMatrix(multiplyMatrix(Solution,B), R_inv);
    Matrix BTP = multiplyMatrix(transposeMatrix(B),Solution);

    ATP = addMatrix(ATP, PA);

    PBR_inv = multiplyMatrix(PBR_inv, BTP);

    ATP = subtraction(ATP, PBR_inv);

    printMatrix(addMatrix(ATP,Q));
    
    deleteMatrix(ATP);
    deleteMatrix(PBR_inv);
    deleteMatrix(Vecs);
    deleteMatrix(Solution);
    deleteMatrix(H_eigen);
    deleteMatrix(H);
    delete[] eigenvalues;
    deleteMatrix(U);
    
    return 0;
}
