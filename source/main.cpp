#include <string>
#include <vector>

#include "data.h"
#include "functions.h"

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

    //Display all the parameters and conditions used for computation
    data->display_parameters();


    /*TODO :
        - Discrétisation Laplacien (space scheme)
        - Time scheme (Euler explicite, Euler implicite, C-N)
        - Résolution de système linéaire AX=b (CG, BiCG, BiCGStab)
        - Ajouter une fonction qui écrit les résultat dans le dossier 'output'
    */


    //Pointer deletion
    delete data;
    return 0;
}