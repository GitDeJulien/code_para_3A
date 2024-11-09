#ifndef _LINEAR_ALGEBRA_H
#define _LINEAR_ALGEBRA_H

#include <string>
#include <iostream>
#include <vector>

#include "data.h"
#include "functions.h"

class LinearAlgebra {

    public:

    LinearAlgebra();

    //Resolve the linear systeme AX=b with LU decomposition and backward, forward methods
    void LU(std::vector<std::vector<double>> A, std::vector<double> X , std::vector<double> b);

    //Resolve the linear systeme AX=b with the Bi-Conjugate Gradient Stabilised methode
    void BiCGStab(std::vector<std::vector<double>> A, std::vector<double> X , std::vector<double> b);

};

#endif