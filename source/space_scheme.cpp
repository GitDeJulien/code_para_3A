#ifndef _SPACE_SCHEME_CPP
#define _SPACE_SCHEME_CPP

#include <iostream>
#include <utility>

#include "space_scheme.h"

SpaceScheme::SpaceScheme(){}

int SpaceScheme::index_MatToVect(Data* data, const int i, const int j)
{
    int Ny(0);
    int l(0);

    Ny = data->Get_Ny();

    l = i*Ny + j;

    return(l);
}

std::pair<int, int> SpaceScheme::index_VectToMat(Data* data, const int l) {

    int Ny(0);
    int i(0);
    int j(0);

    Ny = data->Get_Ny();

    i = l/Ny;
    j = l%Ny;

    return {i,j};
}

std::vector<double> SpaceScheme::Initialize(Data* data, Function* function)
{
    std::vector<double> U0;
    int Nx(0), Ny(0);
    int l(0);
    double x(0.0), y(0.0);

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();

    int N = Nx*Ny;
    U0.resize(N);

    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            l = index_MatToVect(data, i, j);

            x = i*data->Get_hx();
            y = j*data->Get_hy();

            U0[l] = function->InitialCondition(data, x, y);
        }
    }


    return U0;
}

Matrix SpaceScheme::BuildMatrix(Data* data) 
{
    std::pair<int, int> indMat;

    int Nx(0);
    int Ny(0);

    double alpha(0.0);
    double beta(0.0);
    double gamma(0.0);

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();

    int N = Nx*Ny;

    alpha = data->Get_diffusion_coeff()*(2/(data->Get_hx()*data->Get_hx()) + 2/(data->Get_hy()*data->Get_hy()));
    beta = - data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx());
    gamma = - data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy());

    Matrix matrix(N,N);

    if(data->Get_SpaceScheme() == 1){ //Laplacian centered discretisation
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){

                if(i==j){
                    matrix(i,j) = alpha;
                }
                else if((j==i+1 || i==j+1) && i<=Nx){
                    matrix(i,j) = beta;
                }
                else if((j==i+1 || i==j+1) && i>Nx){
                    matrix(i,j) = gamma;
                }
                else {
                    matrix(i,j) = 0.0;
                }
            }
        }
    }
    else {
        std::cerr << "The Space Scheme key" << data->Get_SpaceScheme() << "is not define. Please change it in the 'input/data.dat' file" << std::endl;
        exit(EXIT_FAILURE);
    }
    return matrix;     
} 


std::vector<double> SpaceScheme::SourceTerme(Data* data, Function* function, const double t){

    std::vector<double> S;
    int Nx(0), Ny(0);
    int l;
    double x,y;

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();

    int N = Nx*Ny;
    S.resize(N);

    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            l = index_MatToVect(data, i, j);

            x = i*data->Get_hx();
            y = j*data->Get_hy();

            S[l] = function->SourceFunction(data, x, y, t);
        }
    }

    // TODO!
    // if there is other border condition than homogenus
    // need to modify the source terme in appropriate spotes...
    // make so if, else if, else case do so...

    return S;

}


#endif
