#include "SidewallForceSolver.h"
#include <gsl/gsl_integration.h>
#include <functional>
#include <iostream>
#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <matplotlibcpp.h>
#include "../../math/Integral.h"

namespace plt = matplotlibcpp;
using namespace Hina;
using namespace std;

void SidewallForceSolver::step() 
{

}

float SidewallForceSolver::getDrawbarPull() 
{
    return 0;
}

/*
 * Get the static sinkage of the wheel.
 * @param s: the soil patch.
 * @param wheel: the wheel.
 * @return the static sinkage of the wheel.
 */

float BekkerForceSolver::getStaticSinkage(SoilPatch s, Wheel wheel) 
{
    float theta_s = getStaticContactAngle(s,wheel);
    return wheel.r * ( 1 - cos(theta_s));
}

float SidewallForceSolver::getWheelSinkage() {
    return getStaticSinkage(this->soil,this->wheel) + getDynamicSinkage(this->soil,this->wheel);
}

float BekkerForceSolver::getKineticContactExitAngle(SoilPatch s, Wheel wheel) 
{
    float h_s = getStaticSinkage(s,wheel);
    return acos(1 - (h_s / wheel.r));   
}

float BekkerForceSolver::getKineticContactEntryAngle(SoilPatch s, Wheel wheel) 
{    
    float h_s = getStaticSinkage(s,wheel);
    return -acos(1 - (s.l_s * h_s / wheel.r));
}


/*
 * Compute the side wall resistance force after Schwanghart (1968)
 * @param beta - slip angle
 * @param theta_f - entry angle
 * @return the side wall resistance force in N.
 */
float WheelForceSolver::getSidewallResistance(float beta, float theta_f) 
{
    // passive soil pressure caused by compacting side soil
    double N_q = exp((1.5*M_PI - soil.phi) + tan(soil.phi)) / 2*pow(cos(M_PI/4 + soil.phi/2),2);
    double N_c = (N_q - 1) / tan(soil.phi);
    double N_gamma = 2*(N_q+1)*tan(soil.phi) / (1+0.4*sin(4*soil.phi));
    
    double q = 1; // Fudge parameter (surchage of pressure, e.g. due to bulldozing)
    double sigma_s = soil.gamma*getWheelSinkage()*N_gamma + soil.c*N_c + q*N_q;

    struct integrand_params {
        double h;
        SoilPatch s;
        Wheel w;
    };

    auto integrand = [](double theta, void* params) -> double {
        auto p = *(struct integrand_params*) params;
    
        auto f_x = sqrt(wheel.r - pow(x,2))
        return (p.solver->getTAUX(p.s,p.w, theta) * sin(theta)) + (p.solver->getSigma(p.s,p.w,theta) * cos(theta));
    };

    Integral integral;

    Function fIntegrand(integrand);
    double theta_f = getKineticContactEntryAngle(s, wheel);
    double theta_r = getKineticContactExitAngle(s, wheel);

    auto integrand_param = new integrand_params {
        .solver = this,
        .h = 0,
        .s = s,
        .w = wheel
    };

    IntegralParams params = {
        .integrand = fIntegrand,
        .params = integrand_param,
        .a = theta_f,
        .b = theta_r
    };

    double h = integral.evaluate(params);

    // auto fn = [](double h, void* params) -> double {
    //     auto p = *(struct function_params *) params;
    //     double res = integral.(h, integr_params);
    //     // delete(integr_params);
    //     // return p.wheel.W - (p.wheel.r * p.wheel.b * res);
    //     return -1;
    // };

    delete(integrand_param);

    return -1; 
}