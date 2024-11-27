#ifndef _SPACE_SCHEME_CPP
#define _SPACE_SCHEME_CPP

#include <iostream>
#include <utility>

#include "space_scheme.h"

int SpaceScheme::index_MatToVect(const int i, const int j)
{
    int Nx(0);
    int l(0);

    Nx = data.Get_Nx();

    l = (j-1)*Nx + (i-1);

    return(l);
}

std::pair<int, int> SpaceScheme::index_VectToMat(const int l) {

    int Ny(0);
    int i(0);
    int j(0);

    Ny = data.Get_Ny();

    i = l/Ny+1;
    j = l%Ny+1;

    return {i,j};
}

std::vector<double> SpaceScheme::Initialize(Function* function)
{
    std::vector<double> U0;
    int Nx(0), Ny(0);
    int l(0);
    double x(0.0), y(0.0);

    Nx = data.Get_Nx();
    Ny = data.Get_Ny();

    int N = Nx*Ny;
    U0.resize(N);

    for(int i=1; i<=Nx; ++i){
        for(int j=1; j<=Ny; ++j){
            l = index_MatToVect(i, j);

            x = i*data.Get_hx();
            y = j*data.Get_hy();

            U0[l] = function->InitialCondition(x, y);
        }
    }


    return U0;
}

Matrix SpaceScheme::Initialize2(Function* function)
{
    int Nx = data.Get_Nx();
    int Ny = data.Get_Ny();
    int N = Nx*Ny;
    Matrix U0(Nx,Ny);
    double x, y;

    for(int i=1; i<=Nx; ++i){
        for(int j=1; j<=Ny; ++j){
            x = i*data.Get_hx();
            y = j*data.Get_hy();

            U0(i,j) = function->InitialCondition(x, y);
        }
    }


    return U0;
}

Matrix SpaceScheme::BuildMatrix() 
{
    std::pair<int, int> indMat;

    int Nx(0);
    int Ny(0);

    double alpha(0.0);
    double beta(0.0);
    double gamma(0.0);

    Nx = data.Get_Nx();
    Ny = data.Get_Ny();

    int N = Nx*Ny;

    alpha = data.Get_diffusion_coeff()*(2/(data.Get_hx()*data.Get_hx()) + 2/(data.Get_hy()*data.Get_hy()));
    beta = - data.Get_diffusion_coeff()/(data.Get_hx()*data.Get_hx());
    gamma = - data.Get_diffusion_coeff()/(data.Get_hy()*data.Get_hy());

    Matrix matrix(N,N);

    for (int I=0; I<N; ++I) {
        for (int J=0; J<N; ++J){
            matrix(I,J) = 0.0;
        }
    }

    if(data.Get_SpaceScheme() == 1){ //Laplacian centered discretisation
        for(int i=1; i<=Nx; ++i){
            for(int j=1; j<=Ny; ++j){
                int l = this->index_MatToVect(i, j);

                matrix(l,l) = alpha;
                if (i>1) {
                    matrix(l,l-1) = beta;
                }
                if (i<Nx){
                    matrix(l,l+1) = beta;
                }
                if (j>1) {
                    matrix(l,l-Nx) = gamma;
                }
                if (j<Ny) {
                    matrix(l,l+Nx) = gamma;
                }
            }
        }
    }
    else {
        std::cerr << "The Space Scheme key" << data.Get_SpaceScheme() << "is not define. Please change it in the 'input/data.dat' file" << std::endl;
        exit(EXIT_FAILURE);
    }
    return matrix;     
} 


std::vector<double> SpaceScheme::SourceTerme(Function* function, const double t){

    std::vector<double> S;
    int Nx(0), Ny(0);
    int l;
    double x,y;

    Nx = data.Get_Nx();
    Ny = data.Get_Ny();

    int N = Nx*Ny;
    S.resize(N);

    for(int i=1; i<=Nx; ++i){
        for(int j=1; j<=Ny; ++j){
            l = index_MatToVect(i, j);

            x = i*data.Get_hx();
            y = j*data.Get_hy();

            S[l] = function->SourceFunction(x, y, t);

            if (i==1) { 
                //Derichlet non homogène
                S[l] += data.Get_diffusion_coeff()/(data.Get_hx()*data.Get_hx()) * function->BoundaryCondition_Left(x, y);
            }
            if (i==Nx) { 
                //Derichlet non homogène
                S[l] += data.Get_diffusion_coeff()/(data.Get_hx()*data.Get_hx()) * function->BoundaryCondition_Right(x, y);
            }
            if (j==1) {
                //Derichlet non homogène
                S[l] += data.Get_diffusion_coeff()/(data.Get_hy()*data.Get_hy()) * function->BoundaryCondition_Down(x, y);
            }
            if (j==1) {
                //Derichlet non homogène
                S[l] += data.Get_diffusion_coeff()/(data.Get_hy()*data.Get_hy()) * function->BoundaryCondition_Up(x, y);
            }
        }
    }

    // TODO!
    // if there is other border condition than homogenus
    // need to modify the source terme in appropriate spotes...
    // make so if, else if, else case do so...

    return S;

}

std::pair<Matrix,Matrix> SpaceScheme::BuildMatrix2(){

    int Nx = data.Get_Nx();
    int Ny = data.Get_Ny();
    double hx = data.Get_hx();
    double hy = data.Get_hy();
    double D = data.Get_diffusion_coeff();
    //alpha <-> A, beta <-> B
    double alpha,beta;
    // A matrice gauche, D matrice droite
    Matrix A(Nx,Nx,1),B(Ny,Ny,1);

    alpha = D*(2/(hx*hx));
    beta = D*(2/(hy*hy));

    A = A*alpha;
    B = B*beta;

    if(data.Get_SpaceScheme() == 1){ //Laplacian centered discretisation
        //Construction de A
        for(int i=0; i<Nx; i++){
            for(int j=0; j<Nx; j++){
                if(i == j-1 || i == j+1){
                    A(i,j) = -alpha/2;
                }
            }
        }
        //Construction de B
        for(int i=0; i<Ny; i++){
            for(int j=0; j<Ny; j++){
                if(i == j-1 || i == j+1){
                    B(i,j) = -beta/2;
                }
            }
        }
    }
    else {
        std::cerr << "The Space Scheme key" << data.Get_SpaceScheme() << "is not defined. Please change it in the 'input/data.dat' file" << std::endl;
        exit(EXIT_FAILURE);
    }
    return std::make_pair(A,B);
}

Matrix SpaceScheme::SourceTerm2(Function *function, const double t) {

    int Nx = data.Get_Nx();
    int Ny = data.Get_Ny();
    double hx = data.Get_hx();
    double hy = data.Get_hy();
    double D = data.Get_diffusion_coeff();
    double x,y;

    Matrix S(Nx,Ny);

    for(int i=1; i<=Nx; ++i){
        for(int j=1; j<=Ny; ++j){

            x = i*data.Get_hx();
            y = j*data.Get_hy();

            S(i,j) = function->SourceFunction(x, y, t);

            if (i==1) {
                //Derichlet non homogène
                S(i,j) += D/(hx*hx) * function->BoundaryCondition_Left(x, y);
            }
            if (i==Nx) {
                //Derichlet non homogène
                S(i,j) += D/(hx*hx) * function->BoundaryCondition_Right(x, y);
            }
            if (j==Nx) {
                //Derichlet non homogène
                S(i,j) += D/(hy*hy) * function->BoundaryCondition_Down(x, y);
            }
            if (j==1) {
                //Derichlet non homogène
                S(i,j) += D/(hy*hy) * function->BoundaryCondition_Up(x, y);;
            }
        }
    }

    // TODO!
    // if there is other border condition than homogenus
    // need to modify the source terme in appropriate spotes...
    // make so if, else if, else case do so...

    return S;
}


#endif
