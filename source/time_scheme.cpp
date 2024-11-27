#ifndef _TIME_SCHEME_CPP
#define _TIME_SCHEME_CPP


#include <direct.h>
#include "time_scheme.h"

TimeScheme::TimeScheme(Data data, LinearAlgebra* lin, Function* fct, SpaceScheme* ssch) : data(data), _lin(lin), _fct(fct), _ssch(ssch)
{
    _dt = data.Get_dt();
    _key_TimeScheme = data.Get_TimeScheme();

}

std::vector<double> TimeScheme::EulerExplicite(const Matrix& A, const std::vector<double> Un , const std::vector<double> bn)
{
    if(A.getCols() != static_cast<int>(Un.size()) || Un.size() != bn.size()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    int N = Un.size();

    std::vector<double> Unp1;
    Unp1.resize(N);

    //Identity matrix
    Matrix I = Matrix::Identity(N);

    Unp1 = ((A.ScalarMultiply(-this->_dt)).AddMatrix(I)).MatrixVectorProduct(Un);

    for(int i=0; i<N; ++i) {
        Unp1[i] += this->_dt*bn[i];
    }

    return Unp1;

}

std::vector<double> TimeScheme::EulerImplicite(const Matrix& A, const std::vector<double> Un , const std::vector<double> bnp1)
{

    if(A.getCols() != static_cast<int>(Un.size()) || Un.size() != bnp1.size()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    int N = Un.size();

    std::vector<double> Unp1;
    Unp1.resize(N);

    //Identity matrix
    Matrix I = Matrix::Identity(N);
    Matrix A_star = (A.ScalarMultiply(this->_dt)).AddMatrix(I);
    std::vector<double> U_star;
    U_star.resize(N);

    for(int i=0; i<N; ++i) {
        U_star[i] = Un[i] + this->_dt*bnp1[i];
    }

    //Unp1 = _lin->LU(A_star, U_star);
    Unp1 = _lin->BiCGStab(A_star, U_star, 10000, 1e-6);

    return Unp1;
}

std::vector<double> TimeScheme::CranckNicholson(const Matrix& A, const std::vector<double> Un , const std::vector<double> bn, const std::vector<double> bnp1){

    if(A.getCols() != static_cast<int>(Un.size()) || Un.size() != bn.size() || Un.size() != bnp1.size()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }    

    int N = Un.size();

    std::vector<double> Unp1;
    Unp1.resize(N);

    //Identity matrix
    Matrix I = Matrix::Identity(N);
    Matrix A_star = (A.ScalarMultiply(this->_dt/2.)).AddMatrix(I);
    std::vector<double> U_star;
    U_star.resize(N);

    U_star = ((A.ScalarMultiply(-this->_dt/2.)).AddMatrix(I)).MatrixVectorProduct(Un);

    for(int i=0; i<N; ++i) {
        U_star[i] += this->_dt/2.*(bnp1[i]+bn[i]);
    }

    Unp1 = _lin->BiCGStab(A_star, U_star, 10000, 1e-6);

    return Unp1;
}

std::vector<double> TimeScheme::Advance(const Matrix& A, const std::vector<double> Un, const double tn) 
{
    int N = Un.size();
    std::vector<double> Unp1;
    std::vector<double> Sn;
    std::vector<double> Snp1;
    Unp1.resize(N);
    double tnp1 = tn + data.Get_dt();

    switch (data.Get_TimeScheme())
        {
        case 1:
            Sn.resize(N);
            Sn = _ssch->SourceTerme(_fct, tn);
            Unp1 = this->EulerExplicite(A, Un, Sn);
            break;

        case 2:
            Snp1.resize(N);
            Snp1 = _ssch->SourceTerme(_fct, tnp1);
            Unp1 = this->EulerImplicite(A, Un, Snp1);
            break;
            
        case 3:
            Sn.resize(N);
            Snp1.resize(N);
            Sn = _ssch->SourceTerme(_fct, tn);
            Snp1 = _ssch->SourceTerme(_fct, tnp1);
            Unp1 = this->CranckNicholson(A, Un, Sn, Snp1);
            break;
        
        default:
            std::cout << "The time scheme key " << data.Get_TimeScheme() << " is not valid !" << std::endl;
            exit(EXIT_FAILURE);
    }

    return Unp1;

}

void TimeScheme::SaveSol(const std::vector<double>& sol, const std::string& path, int n) 
{

    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        // Try to create the directory
        if (_mkdir(path.c_str()) == -1) {
            std::cout << "Directory already exist." << std::endl;
        } else {
            std::cout << "Directory " << path << " created successfully." << std::endl;
        }
    }
    std::string n_file = path + "/sol_" + std::to_string(n) + ".vtk";
    std::ofstream solution;
    solution.open(n_file);
    if (solution.is_open()){
        solution << "# vtk DataFile Version 3.0" << std::endl;
        // solution << "#key time scheme: " << data.Get_TimeScheme() << std::endl;
        // solution << "#key space scheme: " << data.Get_SpaceScheme() << std::endl;
        // solution << "#key left/right boundary condition: " << data.Get_key_LeftRightBoundCond() << std::endl;
        // solution << "#key up/down boundary condition: " << data.Get_key_UpDownBoundCond() << std::endl;
        // solution << "#key initial condition: " << data.Get_key_InitialCondition() << std::endl;
        // solution << "#key source terme: " << data.Get_key_SourceTerme() << std::endl;
        solution << "sol" << std::endl;
        solution << "ASCII" << std::endl;
        solution << "DATASET STRUCTURED_POINTS" << std::endl;
        solution << "DIMENSIONS " << data.Get_Nx() << " " << data.Get_Ny() << " " << 1 << std::endl;
        solution << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        solution << "SPACING " << data.Get_hx() << " " << data.Get_hy() << " " << 1 << std::endl;;
        solution << "POINTdata " << data.Get_Nx()*data.Get_Ny() << std::endl;
        solution << "SCALARS sol float" << std::endl;
        solution << "LOOKUP_TABLE default" << std::endl;
        for(int j=0; j<data.Get_Ny(); ++j)
        {
            for(int i=0; i<data.Get_Nx(); ++i)
            {
                solution << sol[i+j*data.Get_Nx()] << " ";
            }
            solution << std::endl;
        }
        solution.close();
    }
    else {
        std::cout << "Error opening results file!" << std::endl;
    }
}

//Fonctions de Nicola

Matrix TimeScheme::Advance(const Matrix& A, const Matrix &B, const Matrix & Un, const double tn)
{
    int rows = Un.getRows();
    int cols = Un.getCols();

    Matrix Unp1 = Un;
    Matrix Sn(rows,cols);
    Matrix Snp1(rows, cols);
    double tnp1 = tn + data.Get_dt();

    switch (data.Get_TimeScheme())
    {
        case 1:
            Sn = _ssch->SourceTerm2(_fct, tn);
            Unp1 = this->EulerExplicit(A, B, Un, Sn);
            break;

        case 2:
            Snp1 = _ssch->SourceTerm2(_fct, tnp1);
            Unp1 = this->EulerImplicit(A, B, Un, Snp1);
            break;

        case 3:
            Sn = _ssch->SourceTerm2(_fct, tn);
            Snp1 = _ssch->SourceTerm2(_fct, tnp1);
            Unp1 = this->CrankNicolson(A, B, Un, Sn, Snp1);
            break;

        default:
            std::cout << "The time scheme key " << data.Get_TimeScheme() << " is not valid !" << std::endl;
            exit(EXIT_FAILURE);
    }

    return Unp1;

}

Matrix TimeScheme::EulerExplicit(const Matrix &A,const Matrix &B, const Matrix & Un, const Matrix & bn) {
    if(A.getCols() != Un.getRows() || Un.getRows() != bn.getRows() || Un.getCols() != bn.getRows()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    int rows = Un.getRows();
    int cols = Un.getCols();

    Matrix Unp1(rows,cols);

    //Identity matrix
    Matrix I(A.getRows(),A.getRows(),1);

    //Returns Unp1
    Unp1 = (A*this->_dt + I)*Un + Un*B*this->_dt + bn*this->_dt;

    return Unp1;
}

Matrix TimeScheme::EulerImplicit(const Matrix &A, const Matrix &B, const Matrix &Un, const Matrix &bnp1) {

    if(A.getCols() != Un.getRows() || Un.getRows() != bnp1.getRows() || Un.getCols() != bnp1.getRows()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    int rows = Un.getRows();
    int cols = Un.getCols();

    Matrix Unp1(rows,cols);

    //Identity matrix
    Matrix I(A.getRows(),A.getRows(),1);

    Matrix A_star = A*this->_dt + I;
    Matrix U_star = Un;

    U_star = Un + bnp1*this->_dt;

    //Unp1 = _lin->LU(A_star, U_star);
    //Unp1 = _lin->BiCGStab(A_star, U_star, 10000, 1e-6);
    Unp1 = _lin->solveLiap(I - A,B,Un + bnp1,10000);

    return Unp1;
}

Matrix TimeScheme::CrankNicolson(const Matrix &A, const Matrix &B,const Matrix &Un, const Matrix &bn, const Matrix &bnp1) {
    if(A.getCols() != Un.getRows() || Un.getRows() != bnp1.getRows() || Un.getCols() != bnp1.getRows() || Un.getRows() != bn.getCols() || Un.getCols() != bn.getRows()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    int rows = Un.getRows();
    int cols = Un.getCols();

    Matrix Unp1(rows,cols);

    //Identity matrix
    Matrix I(A.getRows(),A.getRows(),1);

    Matrix A_star = A*this->_dt + I;
    Matrix U_star = Un;

    U_star = Un + bnp1*this->_dt;

    //Unp1 = _lin->LU(A_star, U_star);
    //Unp1 = _lin->BiCGStab(A_star, U_star, 10000, 1e-6);
    Unp1 = _lin->solveLiap(I - 0.5*A,0.5*B,(I + 0.5*A)*Un + 0.5*Un*B + 0.5*bnp1 + 0.5*bn,10000);

    return Unp1;
}

void TimeScheme::SaveSol(const Matrix & sol, const std::string& path, int n)
{

    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        // Try to create the directory
        if (_mkdir(path.c_str()) == -1) {
            std::cout << "Directory already exist." << std::endl;
        } else {
            std::cout << "Directory " << path << " created successfully." << std::endl;
        }
    }
    std::string n_file = path + "/sol_" + std::to_string(n) + ".vtk";
    std::ofstream solution;
    solution.open(n_file);
    if (solution.is_open()){
        solution << "# vtk DataFile Version 3.0" << std::endl;
        // solution << "#key time scheme: " << data.Get_TimeScheme() << std::endl;
        // solution << "#key space scheme: " << data.Get_SpaceScheme() << std::endl;
        // solution << "#key left/right boundary condition: " << data.Get_key_LeftRightBoundCond() << std::endl;
        // solution << "#key up/down boundary condition: " << data.Get_key_UpDownBoundCond() << std::endl;
        // solution << "#key initial condition: " << data.Get_key_InitialCondition() << std::endl;
        // solution << "#key source terme: " << data.Get_key_SourceTerme() << std::endl;
        solution << "sol" << std::endl;
        solution << "ASCII" << std::endl;
        solution << "DATASET STRUCTURED_POINTS" << std::endl;
        solution << "DIMENSIONS " << data.Get_Nx() << " " << data.Get_Ny() << " " << 1 << std::endl;
        solution << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        solution << "SPACING " << data.Get_hx() << " " << data.Get_hy() << " " << 1 << std::endl;;
        solution << "POINTdata " << data.Get_Nx()*data.Get_Ny() << std::endl;
        solution << "SCALARS sol float" << std::endl;
        solution << "LOOKUP_TABLE default" << std::endl;
        for(int j=0; j<data.Get_Ny(); ++j)
        {
            for(int i=0; i<data.Get_Nx(); ++i)
            {
                solution << sol(i,j)<< " ";
            }
            solution << std::endl;
        }
        solution.close();
    }
    else {
        std::cout << "Error opening results file!" << std::endl;
    }
}

#endif