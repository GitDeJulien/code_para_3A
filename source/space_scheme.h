#ifndef _SPACE_SCHEME_H
#define _SPACE_SCHEME_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include "data.h"
#include "functions.h"
#include "matrix.h"
#include "helpers.h"

class SpaceScheme {

    private:

    double* _sol0; 
    double* _Matrix;
    Data data;

    public:

        //Constructor
        SpaceScheme(Data data): data(data){};

        //Transform matrix index (i,j) in vector index l
        int index_MatToVect(const int i, const int j);

        //Transform matrix index (i,j) in vector index l
        std::pair<int, int> index_VectToMat(const int l);

        //Initialize the solution vector at t=0.0
        std::vector<double> Initialize(Function* function);

        Matrix Initialize2(Function* function);

        //Building the general Matrix M depending on the scheme
        Matrix BuildMatrix();

        std::pair<Matrix,Matrix> BuildMatrix2();

        //Building the Source terme S depending on the border conditions
        std::vector<double> SourceTerme(Function* function, const double t);

        Matrix SourceTerm2(Function* function, const double t);

};





#endif