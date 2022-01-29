#include "WheelForceSolver.h"

using namespace Hina;

void WheelForceSolver::initialize(SoilPatch soil, Wheel wheel)
{ 
    this->soil = soil;
    this->wheel = wheel;
}
