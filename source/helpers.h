//
// Created by UTENTE on 23/11/2024.
//

#define SOURCE_HELPERS_H
#include "matrix.h"
#include <string>
#include <iostream>
#include <vector>

    double operator*(const std::vector<double>& a, const std::vector<double>& b);

    Matrix operator*(const Matrix & A, const Matrix & B);

    Matrix operator*(const Matrix & A, double mu);

    Matrix operator*(double mu,const Matrix & A);

    std::vector<double> operator*(const Matrix & A, const std::vector<double> & b);

    double norm(const std::vector<double> & a);

    double norm(const Matrix & A);

    Matrix operator+(const Matrix& A, const Matrix& B);

    Matrix operator-(const Matrix& A, const Matrix& B);

    Matrix trMatrix(const Matrix & A);

    std::vector<double> computeAdiCoefficients(int n, double diagA, double subdiagA, double diagB, double subdiagB, int numCoefficients);

    void computeSpectralBounds(int n, double diag, double subdiag, double &lambda_min, double &lambda_max);


