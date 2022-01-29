#ifndef WHEEL_FORCE_SOLVER
#define WHEEL_FORCE_SOLVER

#include "../body/Wheel.h"
#include "../environment/soil/SoilPatch.h"

namespace Hina {
    
    class WheelForceSolver {

    private:
        Wheel wheel;
        SoilPatch soil;
        
    public:
        WheelForceSolver() {};
        ~WheelForceSolver() {};

        void initialize(SoilPatch soil, Wheel wheel);
        virtual void step() {};
        virtual float getDrawbarPull();
        virtual float getWheelSinkage();
    };
}

#endif