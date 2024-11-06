#include <string>
#include <vector>

#include "data.h"
#include "functions.h"

int main(int argc, char** argv) {
    std::string filename = "../input/data.dat";

    //Pointer to Data class
    Data* data = new Data(filename);

    //Pointer to Function class
    Function* function = new Function(data);

    //Display all the parameters and conditions used for computation
    data->display_parameters();

    //Test function InitialCondition
    std::cout << "Test InitialCondition (test result is 4): " << function->InitialCondition(data, 2.0, 2.0, 3.0) << std::endl;


    /*TODO :
        - Makefile
        - Discrétisation Laplacien (space scheme)
        - Time scheme (Euler explicite, Euler implicite, C-N)
        - Résolution de système linéaire AX=b (CG, BiCG, BiCGStab)
        - Ajouter une fonction qui écrit les résultat dans le dossier 'output'
    */


    //Pointer deletion
    delete data;
    return 0;
}