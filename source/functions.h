#ifndef _FUNCTION_H

#include <fstream>
#include <iostream>
#include <cmath>

#include "data.h"

class Function {

    // private:


    //     double _Lx;
    //     double _Ly;
    //     int _key_InitialCondition;
    //     int _key_SourceTerme;
    //     int _key_LeftRightBoundCond;
    //     int _key_UpDownBoundCond;



    public:
        //Constructor
        Function(Data* data);

        // Initial condition
        const double InitialCondition(const Data* data, const double x, const double y, const double t) const;

        // Source terme
        const double SourceFunction(const Data* data, const double x, const double y, const double t) const;

        // Exacte solution if it known (uusefull for validation)
        const double ExactSolution(const Data* data, const double x, const double y, const double t) const;

        //Boundary condition
        const double BoundaryCondition_g(const Data* data, const double x, const double y, const double t) const;

        const double BoundaryCondition_h(const Data* data, const double x, const double y, const double t) const;

        


};

#define _FUNCTION_H
#endif