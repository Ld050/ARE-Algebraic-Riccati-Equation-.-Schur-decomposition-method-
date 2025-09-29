#include "matrix.hpp"
#include "gauss.hpp"
#include <iostream>
//A: n x n
//B: n x m
//Q: n x n
//R: m x m

//Построение гамильтоновой структуры.
Matrix constructHamiltonian(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R) {
    int n = A.rows;
    if((A.rows != B.rows) || (A.rows != Q.rows) || (Q.rows != Q.cols) || (B.cols != R.rows) || (R.rows != R.cols)){
        std::cout << "Unable to construct Hamiltonian. Check matrix dimensions.\n" << std::endl;
        std::cout << "A dimension\n" << "A: " << A.rows << " x " << A.cols << std:: endl;
        std::cout << "B dimension:\n" << "B: " << B.rows << " x " << B.cols << std:: endl;
        std::cout << "Q dimension:\n" << "Q: " << Q.rows << " x " << Q.cols << std:: endl;
        std::cout << "R dimension:\n" << "R: " << R.rows << " x " << R.cols << std:: endl;
        exit(1);
    }
    Matrix H = createMatrix(2 * n, 2 * n);
    Matrix R_cpy = createMatrix(R.rows, R.cols);
    for(int i = 0; i < R.rows; i++){
        for(int j = 0; j < R.cols; j++){
            R_cpy.data[i][j] = R.data[i][j];
        }
    }
    // A и -A^T 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H.data[i][j] = A.data[i][j];
            H.data[n + i][n + j] = -A.data[j][i]; // -A^T 
        }
    }

    // -BR^{-1}B^T 
    Matrix R_inv = invertMatrix(R_cpy);
    Matrix BR_inv = multiplyMatrix(B, R_inv);
    Matrix BR_inv_B_t = multiplyMatrix(BR_inv, transposeMatrix(B));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H.data[i][n + j] = -BR_inv_B_t.data[i][j];
        }
    }

    // -Q 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H.data[n + i][j] = -Q.data[i][j];
        }
    }

    
    deleteMatrix(BR_inv);
    deleteMatrix(BR_inv_B_t);

    return H;
}
