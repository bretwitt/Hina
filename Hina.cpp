#include <iostream>
#include "src/core/solver/BekkerForceSolver/BekkerForceSolver.h"
#include <math.h>
#include <vector>
#include <matplotlibcpp.h>

using namespace Hina;
namespace plt = matplotlibcpp;


#define DEG_TO_RAD 0.0174533

int main() {

    BekkerForceSolver sol;
    SoilPatch s;
    Wheel w;
    w.r = 0.09;
    w.b = 0.11;
    w.W = 64.68f;
    w.beta = 5.0f;
    w.s = 0.5f;

    sol.initialize(s,w);
    
    s.init_parameters(0.8f*1000, 37.2f * DEG_TO_RAD, 1370, 814000, 1, 0.40f, 0.15f, 1600, 0.9f, (0.043 * 5.0f * DEG_TO_RAD) + 0.036f, (0.020f * 5.0f * DEG_TO_RAD) + 0.013f);


    float contact_angles = sol.getStaticContactAngle(s,w);
    std::cout << contact_angles << std::endl;
    
    float sinkage = sol.getStaticSinkage(s,w);
    std::cout << sinkage << std::endl;

    float dyn_sinkage = sol.getDynamicSinkage(s,w);
    std::cout << dyn_sinkage << std::endl;

    // std::vector<double> x,y;

    // for(double i = -M_PI_2; i < M_PI_2; i+=0.01) {
    //     double _y = sol.getTAUX(s, w, i, 0.016);
    //     x.push_back(i);
    //     y.push_back(_y);
    //     std::cout << _y << std::endl;
    // }    
    // plt::plot(x, y);
    // plt::title("Tau");
    // plt::show();


}