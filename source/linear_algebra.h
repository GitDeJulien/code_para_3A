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

    LinearAlgebra();

    //Resolve the linear systeme AX=b with LU decomposition and backward, forward methods
    std::vector<double> LU(const Matrix& A, const std::vector<double> b);

    //Resolve the linear systeme AX=b with the Bi-Conjugate Gradient Stabilised methode
    std::vector<double> BiCGStab(const Matrix& A, const std::vector<double> b);

};

#endif