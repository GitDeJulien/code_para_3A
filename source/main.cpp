#include <string>
#include <vector>
#include<time.h>
//#include <mscat.h>

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
    Data data(filename);

    //Pointer to Function class
    Function function(data);

    //Pointer to Linear Algebra class
    LinearAlgebra lin = LinearAlgebra();

    //Pointer to Space Scheme class
    SpaceScheme ssch = SpaceScheme(data);

    //Pointer to Time Scheme class
    TimeScheme tsch = TimeScheme(data, &lin, &function, &ssch);

    //Pointer to Matrix class
    //Matrix* matrix = new Matrix();

    //Display all the parameters and conditions used for computation
    data.display_parameters();

    int Nx(0);
    int Ny(0);
    int N(0);

    Nx = data.Get_Nx();
    Ny = data.Get_Ny();
    N = Nx*Ny;

    Matrix A(N,N);
    Matrix Un(Nx,Ny);
    Matrix Sn(Nx,Ny);

    //Laplacian matrix discretisation
    std::pair<Matrix,Matrix> matrices = ssch.BuildMatrix2();

    //Initial solution
    Un = ssch.Initialize2(&function);

    Matrix Unp1 = Un;

    tsch.SaveSol(Un, data.Get_outputPath(), 0);


    double tn = data.Get_t0();
    double nb_iteration = data.Get_niter();
    double dt = data.Get_dt();

    for(int iter = 1; iter<nb_iteration; ++iter) {

        //Advance of a time step with the chosen time scheme
        Unp1 = tsch.Advance(matrices.first, matrices.second, Un, tn);

        //Download result in vtk files
        tsch.SaveSol(Unp1, data.Get_outputPath(), iter);

        //Update
        Un = Unp1;
        tn += dt;
        std::cout << "tn: " << tn << std::endl;
    }

    Matrix U_exact(Nx,Ny);
    double x;
    double y;
    for(int i=1; i<=Nx; ++i){
        for(int j=1; j<=Ny; ++j){
            x = i*data.Get_hx();
            y = j*data.Get_hy();
            U_exact(i-1,j-1) = function.ExactSolution(x, y);
        }
    }

    tsch.SaveSol(U_exact, "output/Exact2", 0);

    //TODO :

        //- Commencer à réfléchir à une stratégie de parallélisation

    return 0;
}


//############################//
//########## TEST ############//
//############################//
/*
// Test Liapunov
int main(){
    Matrix A(3,3,1);
    Matrix B(3,3,1);
    Matrix C(3,3,1);
    B = 10*B;
    A = 2*A;
    C(0,1) = -1;
    C(1,0) = 3;
    B(0,1) = -1;
    B(1,0) = -1;
    A(1,0) = -1;
    A(0,1) = -1;
    A(1,2) = -1;
    A(2,1) = -1;
    B(2,1) = -1;
    B(1,2) = -1;
    A.print();
    B.print();
    C.print();
    LinearAlgebra lin;
    Matrix X = lin.solveLiap(A,B,C,1000);
    X.print();
}*/
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
