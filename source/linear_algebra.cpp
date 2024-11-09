#ifndef _LINEAR_ALGEBRA_CPP
#define _LINEAR_ALGEBRA_CPP

#include "linear_algebra.h";


LinearAlgebra::LinearAlgebra(){};

void LinearAlgebra::LU(std::vector<std::vector<double>> A, std::vector<double> X , std::vector<double> b){

    int N = b.size();

    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;

    U = A;

    if (A.size() != N) {
        std::cerr << "Error: dimensions don't agree!" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
        L[i][i] = 1.;
    }

    //building L and U via gauss elimination
    for (int k = 0; k < N; ++k) {
        for (int i = k+1; i < N; i++) {
            double multiplier = U[i][k]/U[k][k];
            L[i][k] = multiplier;
            U[i][k] = 0.;
            for (int j = k+1; j < N; j++) {
                U[i][j] -= multiplier*U[k][j];
            }
        }
    }

    
    if (U.size() != N) {
        std::cerr << "Error: dimensions don't agree!" << std::endl;
        exit(EXIT_FAILURE);
    }

    //backward solving Ux=y
    for (int i = 1; i <= N; ++i) {
        double s = 0.;
        for (int j = 1; j < i; ++j){
            s += U[N-i][N-j]*X[N-j];
        }
        X[N-i] = (b[N-i] - s)/U[N-i][N-i];
    }

    //forward solving Ly = b
    for (int i = 0; i < N; ++i) {
        double s = 0.;
        for (int j = 0; j < i; ++j)
            s += L[i][j]*X[j];
        X[i] = (b[i] - s)/L[i][i];
    }

}

void LinearAlgebra::BiCGStab(std::vector<std::vector<double>> A, std::vector<double> X , std::vector<double> b){
    
}


#endif