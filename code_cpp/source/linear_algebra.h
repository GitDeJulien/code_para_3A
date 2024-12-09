#ifndef _LINEAR_ALGEBRA_H
#define _LINEAR_ALGEBRA_H

#include <string>
#include <iostream>
#include <vector>

#include "data.h"
#include "functions.h"
#include "matrix.h"

class LinearAlgebra {

public:

    double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
            result += a[i] * b[i];
        return result;
    }

    double norm(const std::vector<double>& v) {
        return std::sqrt(dot(v, v));
    }

    //Operator overload
    Matrix& operator*(const Matrix & b) const{
        int n = this->operator*()
    };

    LinearAlgebra();

    //Resolve the linear systeme AX=b with LU decomposition and backward, forward methods
    std::vector<double> LU(const Matrix& A, const std::vector<double> b);

    //Resolve the linear systeme AX=b with the Bi-Conjugate Gradient Stabilised methode
    std::vector<double> BiCGStab(const Matrix& A, const std::vector<double> b, int maxIterations, double tol);

};

#endif