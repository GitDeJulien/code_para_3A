#ifndef _DATA_H
#define _DATA_H

#include <string>
#include <vector>
#include <iostream>

class Data {

    private:
        
        //Space
        double _Lx;
        double _Ly;
        double _hx;
        double _hy;

        //Time
        double _t0;
        double _tfinal;
        double _dt;

        //Time scheme
        int _key_TimeScheme;

        //Space scheme
        int _key_SpaceScheme;

        //Boundary conditions key
        int _key_LeftRightBoundCond;
        int _key_UpDownBoundCond;

        //Source terme key
        int _key_SourceTerme;

        //Initial condition key
        int _key_InitialCondition;

        //Other
        double _diffusionCoeff;

        


    public:

        // Constructeur
        Data(std::string file_name);

        
        const double & Get_diffusion_coeff() const {return _diffusionCoeff;};
        const double & Get_Lx() const {return _Lx;};
        const double & Get_Ly() const {return _Ly;};
        const double & Get_hx() const {return _hx;};
        const double & Get_hy() const {return _hy;};
        const double & Get_t0() const {return _t0;};
        const double & Get_tfinal() const {return _tfinal;};
        const double & Get_dt() const {return _dt;};
        const int & Get_TimeScheme() const {return _key_TimeScheme;};
        const int & Get_SpaceScheme() const {return _key_SpaceScheme;};
        const int & Get_key_LeftRightBoundCond() const {return _key_LeftRightBoundCond;};
        const int & Get_key_UpDownBoundCond() const {return _key_UpDownBoundCond;};
        const int & Get_key_SourceTerme() const {return _key_SourceTerme;};
        const int & Get_key_InitialCondition() const {return _key_InitialCondition;};


        void display_parameters() const;


};

#endif
