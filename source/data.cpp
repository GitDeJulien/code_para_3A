#ifndef _DATA_CPP

#include "data.h"
#include <fstream>
#include <iostream>
//#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <string>



// Constructeur
Data::Data(std::string file_name)
{

    //Data file reading
    std::ifstream file(file_name);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << file_name << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Ignore empty lines and labels (lines starting with '#')
        if (line.empty() || line[0] == '#') continue;

        // Parse key-value pairs
        size_t equal_pos = line.find('=');
        if (equal_pos != std::string::npos) {
            std::string key = line.substr(0, equal_pos - 1);
            std::string string_value = line.substr(equal_pos + 1);
            //double double_value = std::stod(line.substr(equal_pos + 1));
            double double_value = atof(line.substr(equal_pos + 1).c_str());
            int int_value = atoi(line.substr(equal_pos + 1).c_str());


            // Assign values based on the key
            if (key == "Lx") this->_Lx = double_value;
            else if (key == "Ly") this->_Ly = double_value;
            else if (key == "hx") this->_hx = double_value;
            else if (key == "hy") this->_hy = double_value;
            else if (key == "t0") this->_t0 = double_value;
            else if (key == "tfinal") this->_tfinal = double_value;
            else if (key == "dt") this->_dt = double_value;
            else if (key == "D") this->_diffusionCoeff = double_value;
            else if (key == "key_time_scheme") this->_key_TimeScheme = int_value;
            else if (key == "key_space_scheme") this->_key_SpaceScheme = int_value;
            else if (key == "key_LeftRightBoundCond") this->_key_LeftRightBoundCond = int_value;
            else if (key == "key_UpDownBoundCond") this->_key_UpDownBoundCond = int_value;
            else if (key == "key_SourceTerme") this->_key_SourceTerme = int_value;
            else if (key == "key_InitialCondition") this->_key_InitialCondition = int_value;
            else std::cerr << "Error: This parameters is not valid : " << key << std::endl;

        }
    }

    file.close();

}

const void Data::display_parameters() const {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "-----------------PARAMETERS------------------" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Space Parameters:" << std::endl;
    std::cout << "Lx = " << _Lx << ", Ly = " << _Ly << std::endl;
    std::cout << "hx = " << _hx << ", hy = " << _hy << std::endl;
    std::cout << "---------------------------------------------" << std::endl;


    std::cout << "Time Parameters:" << std::endl;
    std::cout << "t0 = " << _t0 << ", tfinal = " << _tfinal << ", dt = " << _dt << std::endl;
    std::cout << "---------------------------------------------" << std::endl;



    std::cout << "Boundary condition:" << std::endl;
    if (_key_LeftRightBoundCond == 1) std::cout << "h = 0 (left/right)" << std::endl;
    else if (_key_LeftRightBoundCond == 2) std::cout << "h = 1 (left/right)" << std::endl;
    else if (_key_LeftRightBoundCond == 3) std::cout << "h = sin(x)+cos(y) (left/right)" << std::endl;
    else std::cerr << "Error: The LeftRight boundary condition key" << _key_LeftRightBoundCond << "is not referenced" << std::endl;

    if (_key_UpDownBoundCond == 1) std::cout << "g = 0 (up/down)" << std::endl;
    else if (_key_UpDownBoundCond == 2) std::cout << "g = 1 (up/down)" << std::endl;
    else if (_key_UpDownBoundCond == 3) std::cout << "g = sin(x)+cos(y) (up/down)" << std::endl;
    else std::cerr << "Error: The UpDown boundary condition key" << _key_UpDownBoundCond << "is not referenced" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    

    if (_key_TimeScheme == 1) std::cout << "Time scheme: " << "Explicite Euler" << std::endl;
    else if (_key_TimeScheme == 2) std::cout << "Time scheme: " << "Implicte Euler " << std::endl;
    else if (_key_TimeScheme == 3) std::cout << "Time scheme: " << "Crank-Nicholson" << std::endl;
    else std::cerr << "Error: The time scheme key" << _key_TimeScheme << "is not referenced" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    if (_key_SpaceScheme == 1) std::cout << "Space scheme: " << "Laplacian Centered" << std::endl;
    else std::cerr << "Error: The time scheme key" << _key_SpaceScheme << "is not referenced" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    if (_key_SourceTerme == 1) std::cout << "Source terme: " << "f = 2*(x-x*x+y-y*y)" << std::endl;
    else if (_key_SourceTerme == 2) std::cout << "Source terme: " << "f = sin(x) + cos(y)" << std::endl;
    else if (_key_SourceTerme == 3) std::cout << "Source terme: " << "f = exp(-(x-Lx/2)^2)*exp(-(y-_Ly/2)^2)*cos(pi/2*t)" << std::endl;
    else std::cerr << "Error: The source terme key" << _key_SourceTerme << "is not referenced" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    if (_key_InitialCondition == 1) std::cout << "Initial condition: " << "u0 = ExactSol(x,y,0)" << std::endl;
    else if (_key_InitialCondition == 2) std::cout << "Initial condition: " << "u0 = sin(x) + cos(y)" << std::endl;
    else if (_key_InitialCondition == 3) std::cout << "Initial condition: " << "u0 = ..." << std::endl;
    else std::cerr << "Error: The initial condition key" << _key_InitialCondition << "is not referenced" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;


    std::cout << "Diffusion coefficient:" << std::endl;
    std::cout << "D = " << _diffusionCoeff << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout <<  std::endl;

}


#define _DATA_CPP
#endif