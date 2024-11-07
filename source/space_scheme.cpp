#ifndef _SPACE_SCHEME_CPP
#define _SPACE_SCHEME_CPP

#include <iostream>
#include <utility>

#include "space_scheme.h"

SpaceScheme::SpaceScheme(){}

int SpaceScheme::index_MatToVect(Data* data, int i, int j)
{
    int Ny(0);
    int l(0);

    Ny = (int)floor(data->Get_Ly()/data->Get_hy()) + 1;
    l = i*Ny + j;

    return(l);
}

std::pair<int, int> SpaceScheme::index_VectToMat(Data* data, int l) {

    int Ny(0);
    int i(0);
    int j(0);

    Ny = (int)floor(data->Get_Ly()/data->Get_hy()) + 1;
    i = l/Ny;
    j = l%Ny;

    return {i,j};
}

std::vector<double> SpaceScheme::Initialize(Data* data, Function* function)
{
    
}

std::vector<std::vector<double>> SpaceScheme::BuildMatrix(Data* data) 
{
    int i(0);
    int j(0);
    std::pair<int, int> indMat;

    int Nx(0);
    int Ny(0);

    double alpha(0.0);
    double beta(0.0);
    double gamma(0.0);

    Nx = (int)floor(data->Get_Lx()/data->Get_hx()) + 1;
    Ny = (int)floor(data->Get_Ly()/data->Get_hy()) + 1;

    alpha = data->Get_diffusion_coeff()*(2/data->Get_hx() + 2/data->Get_hy());
    beta = - data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx());
    gamma = - data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy());

    std::vector<std::vector<double>> matrix;

    for(int l=0; l<Nx*Ny; ++l) {

        indMat = index_VectToMat(data, l);

        //Remplir matrix
    } 



}




























#endif
