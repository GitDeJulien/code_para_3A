//
// Created by UTENTE on 13/11/2024.
//

#ifndef SOURCE_HELPERS_H
#define SOURCE_HELPERS_H

#include <string>
#include <iostream>
#include <vector>
#include "matrix.h"
#include "helpers.h"
#include "math.h"


double operator*(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
    return result;
}

Matrix operator*(const Matrix & A, const Matrix & B) {
    if (A.getCols() == B.getRows()){
        Matrix C(A.getRows(),B.getCols());
        for (size_t j = 0; j < A.getRows(); j++){
            for (size_t i = 0; i < B.getCols(); ++i) {
                for (size_t k = 0; k < A.getCols(); k++){
                    C(i,j) += A(i,k)*B(k,j);
                }
            }
        }
        return C;
    }
}

Matrix operator*(const Matrix & A, double mu) {
    Matrix B(A.getRows(),A.getCols());
    for (size_t j = 0; j < A.getRows(); j++){
        for (size_t i = 0; i < B.getCols(); ++i) {
            B(i,j) = mu*A(i,j);
        }
    }
    return B;
}

Matrix operator*(double mu,const Matrix & A) {
    Matrix B(A.getRows(),A.getCols());
    for (size_t j = 0; j < A.getRows(); j++){
        for (size_t i = 0; i < B.getCols(); ++i) {
            B(i,j) = mu*A(i,j);
        }
    }
    return B;
}

std::vector<double> operator*(const Matrix & A, const std::vector<double> & b) {
    if (A.getCols() == b.size()){
        std::vector<double> c;
        c.resize(b.size());
        for (size_t j = 0; j < A.getRows(); j++){
            for (size_t i = 0; i < b.size(); ++i) {
                c[j] += A(j,i)*b[i];
            }
        }
        return c;
    }
}

double norm(const std::vector<double> & a){
    double n = 0;
    for(size_t i = 0; i < a.size(); i++){
        n += a[i]*a[i];
    }
    return n;
}

double norm(const Matrix & A){
    double n = 0;
    for(size_t i = 0; i < A.getRows(); i++){
        for(size_t j=0; j< A.getCols(); j++){
            n += A(i,j)*A(i,j);
        }
    }
    return sqrt(n);
}

Matrix operator+(const Matrix& A, const Matrix& B){
    if(A.getCols()==B.getCols() && A.getRows()==B.getRows()){
        Matrix C(A.getRows(),A.getCols());
        for(int i = 0; i<A.getRows(); i++){
            for(int j = 0; j- A.getCols(); j++){
                C(i,j) = A(i,j) + B(i,j);
            }
        }
        return C;
    }
}

Matrix operator-(const Matrix& A, const Matrix& B){
    if(A.getCols()==B.getCols() && A.getRows()==B.getRows()){
        Matrix C(A.getRows(),A.getCols());
        for(int i = 0; i<A.getRows(); i++){
            for(int j = 0; j- A.getCols(); j++){
                C(i,j) = A(i,j) - B(i,j);
            }
        }
        return C;
    }
}

Matrix trMatrix(const Matrix & A){
    int r = A.getRows();
    int c = A.getCols();
    Matrix B(c,r);
    for(size_t i=0; i<A.getRows();i++){
        for(size_t j=0; j<A.getCols();j++){
            B(i,j) = A(j,i);
        }
    }
    return B;
}

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// Funzione per calcolare gli autovalori estremi di una matrice tridiagonale
void computeSpectralBounds(int n, double diag, double subdiag, double &lambda_min, double &lambda_max) {
    // Gli autovalori di una matrice tridiagonale simmetrica sono noti:
    // lambda_k = diag + 2 * subdiag * cos(pi * k / (n + 1)) per k = 1, ..., n
    lambda_min = diag + 2 * subdiag * cos(M_PI * n / (n + 1)); // Estremo minimo
    lambda_max = diag + 2 * subdiag * cos(M_PI * 1 / (n + 1)); // Estremo massimo
}

// Funzione per calcolare i coefficienti ottimali ADI
std::vector<double> computeAdiCoefficients(int n, double diagA, double subdiagA, double diagB, double subdiagB, int numCoefficients) {
    double lambdaA_min, lambdaA_max, lambdaB_min, lambdaB_max;

    // Calcola i limiti spettrali delle matrici A e B
    computeSpectralBounds(n, diagA, subdiagA, lambdaA_min, lambdaA_max);
    computeSpectralBounds(n, diagB, subdiagB, lambdaB_min, lambdaB_max);

    // Combina gli intervalli spettrali
    double alpha = std::min(lambdaA_min, lambdaB_min);
    double beta = std::max(lambdaA_max, lambdaB_max);

    // Calcolo dei coefficienti ottimali usando gli zeri dei polinomi di Chebyshev
    std::vector<double> coefficients(numCoefficients);
    for (int k = 0; k < numCoefficients; ++k) {
        double theta = (2 * k + 1) * M_PI / (2.0 * numCoefficients);
        coefficients[k] = alpha + (beta - alpha) * (1 + cos(theta)) / 2.0;
    }

    return coefficients;
}



#endif //SOURCE_HELPERS_H
