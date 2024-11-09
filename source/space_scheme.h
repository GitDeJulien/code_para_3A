#ifndef _SPACE_SCHEME_H
#define _SPACE_SCHEME_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include "data.h"
#include "functions.h"

class SpaceScheme {

    private:

    std::vector<double> _sol0; 
    std::vector<std::vector<double>> _Matrix;

    public:

        //Constructor
        SpaceScheme();

        //Transforme matrix index (i,j) in vector index l
        int index_MatToVect(Data* data, int i, int j);

        //Transforme matrix index (i,j) in vector index l
        std::pair<int, int> index_VectToMat(Data* data, int l);

        //Initialize the solution vector at t=0.0
        std::vector<double> Initialize(Data* data, Function* function);

        //Building the general Matrix M depending on the scheme
        std::vector<std::vector<double>> BuildMatrix(Data* data);

        //Building the Source terme S depending on the border conditions
        std::vector<double> SourceTerme(Data* data, Function* function);

};





#endif