#ifndef _LINEAR_ALGEBRA_H
#define _LINEAR_ALGEBRA_H

#include <string>
#include <iostream>
#include <vector>

#include "data.h"
#include "functions.h"
#include "matrix.h"
#include "helpers.h"

class LinearAlgebra {

public:

    double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
            result += a[i] * b[i];
        return result;
    }

    /*
    double norm(const std::vector<double>& v) {
        return std::sqrt(dot(v, v));
    }
     */

    LinearAlgebra();

    //Resolve the linear systeme AX=b with LU decomposition and backward, forward methods
    std::vector<double> LU(const Matrix& A, const std::vector<double> b);

    //Solve through LU decomposition the system AX=B
    Matrix LU(const Matrix& A, const Matrix& B);

    //Resolve the linear systeme AX=b with the Bi-Conjugate Gradient Stabilised methode
    std::vector<double> BiCGStab(const Matrix& A, const std::vector<double> b, int maxIterations, double tol);

    Matrix BiCGStab(const Matrix& A, const Matrix & B, int maxIterations, double tol);

    //Solves Liapunov problem: Find X s.t. AX - BX = C
    Matrix solveLiap(const Matrix & A, const Matrix & B, const Matrix & C, int nit);

};

#endif