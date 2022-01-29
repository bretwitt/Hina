#include <iostream>
#include "src/core/solver/BasicForceSolver/BasicForceSolver.h"

using namespace Hina;

int main() {

    BasicForceSolver sol;
    SoilPatch s;
    Wheel w;
    w.r = 1;
    w.b = 3;

    s.init_parameters(0.8f*1000, 37.2f, 1370, 814000, 1, 0.40f, 0.15f, 1600, 1, (0.043 * 0.0f) + 0.036f, (0.020f * 0.0f) + 0.013f);

    float contact_angles = sol.getStaticContactAngle(s,w);
    std::cout << contact_angles << std::endl;


}