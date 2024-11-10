#include <string>
#include <vector>
#include<time.h>

#include "data.h"
#include "functions.h"
#include "matrix.h"
#include "linear_algebra.h"
#include "space_scheme.h"
#include "time_scheme.h"



int main(int argc, char** argv) {

    if (argc < 2)
    {
        std::cout << "Please, enter the name of your data file." << std::endl;
        exit (EXIT_FAILURE);
    }

    //std::string filename = "./input/data.dat";
    std::string filename = argv[1];

    //Pointer to Data class
    Data* data = new Data(filename);

    //Pointer to Function class
    Function* function = new Function();

    //Pointer to Linear Algebra class
    LinearAlgebra* lin = new LinearAlgebra();

    //Pointer to Space Scheme class
    SpaceScheme* ssch = new SpaceScheme();

    //Pointer to Time Scheme class
    TimeScheme* tsch = new TimeScheme(data, lin, function, ssch);

    //Pointer to Matrix class
    //Matrix* matrix = new Matrix();

    //Display all the parameters and conditions used for computation
    data->display_parameters();

    int Nx(0);
    int Ny(0);
    int N(0);

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();
    N = Nx*Ny;

    Matrix A(N,N);
    std::vector<double> Un;
    std::vector<double> Unp1;
    std::vector<double> Sn;

    //Laplacian matrix discretisation
    A = ssch->BuildMatrix(data);

    //Initial solution
    Un = ssch->Initialize(data, function);

    //Source terme
    Sn = ssch->SourceTerme(data, function, 0.0);

    Unp1 = Un;

    double tn = data->Get_t0();
    double nb_iteration = data->Get_niter();
    double dt = data->Get_dt();

    for(int iter = 0; iter<nb_iteration; ++iter) {

        //Advance of a time step with the chosen time scheme
        Unp1 = tsch->Advance(A, Un, tn);

        //Download result in vtk files
        tsch->SaveSol(Unp1, "output/cas1", iter);

        //Update
        Un = Unp1;
        tn += dt;
        std::cout << "tn: " << tn << std::endl;
    }



    /*TODO :

        - Coder le BiCGStab
        - Comprendre pourquoi ça ne marche pas avec Euler Implicite
        - Faire une sortie de la solution exacte pour pouvoir comparer en coupe
        - Commencer à réfléchir à une stratégie de parallélisation
    */


    //Pointer deletion
    delete data, delete function, delete lin, delete ssch, delete tsch;
    return 0;
}


//############################//
//########## TEST ############//
//############################//

/*
TEST MATRIX CLASS

Matrix A(2, 3); // 2x3 matrix
A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;

std::vector<double> x = {1, 2, 3};
std::vector<double> result_vector = A.MatrixVectorProduct(x);

std::cout << "Matrix-Vector Product:" << std::endl;
for (double val : result_vector) {
    std::cout << val << " ";
}
std::cout << std::endl;

Matrix B(3, 2); // 3x2 matrix
B(0, 0) = 7; B(0, 1) = 8;
B(1, 0) = 9; B(1, 1) = 10;
B(2, 0) = 11; B(2, 1) = 12;

Matrix C = A.MatrixMatrixProduct(B);

std::cout << "Matrix-Matrix Product:" << std::endl;
for (int i = 0; i < C.rows; ++i) {
    for (int j = 0; j < C.cols; ++j) {
        std::cout << C(i, j) << " ";
    }
    std::cout << std::endl;
}

Matrix A(2, 2); // 2x2 matrix
A(0, 0) = 1; A(0, 1) = 2;
A(1, 0) = 3; A(1, 1) = 4;

Matrix B(2, 2); // Another 2x2 matrix
B(0, 0) = 5; B(0, 1) = 6;
B(1, 0) = 7; B(1, 1) = 8;

Matrix C = A.Add(B);

std::cout << "Matrix Sum:" << std::endl;
for (int i = 0; i < C.rows; ++i) {
    for (int j = 0; j < C.cols; ++j) {
        std::cout << C(i, j) << " ";
    }
    std::cout << std::endl;
}

*/
