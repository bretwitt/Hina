#include "BasicForceSolver.h"
#include <gsl/gsl_integration.h>
#include <functional>
#include <iostream>

using namespace Hina;
using namespace std;
using namespace std::placeholders; 

void BasicForceSolver::step() 
{

}

float BasicForceSolver::getDrawbarPull() 
{
    return 0;
}   

float BasicForceSolver::getWheelSinkage() 
{
    return 0;
}

float BasicForceSolver::getStaticSinkage(SoilPatch s, Wheel wheel) 
{
    return 0;
}

float BasicForceSolver::getStaticContactAngle(SoilPatch s, Wheel wheel) 
{
    float theta_s = -1;
    float k = pow(wheel.r, s.n + 1) * (s.k_c + (s.k_phi * wheel.b));

    auto integrand = [](SoilPatch s, float theta, float theta_s) -> double {
        return cos(theta) - (pow(cos(theta_s), s.n) * cos(theta)); 
    };

    using integrand_ptr = double(*)(SoilPatch,float,float);
    integrand_ptr integr = integrand;

    auto integral = [] (SoilPatch s, float theta_s, integrand_ptr integrand) {
        // Integrate over theta from 0 to pi
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        double result, error;
        gsl_function F;
        F.function = integrand;
        F.params = &s;
        gsl_integration_qags (&F, 0, M_PI, 0, 1e-7, 1000, w, &result, &error);
        gsl_integration_workspace_free (w);
        return result;
    };
    return integral(s,0,integr);
}

float BasicForceSolver::getSigma(SoilPatch s, Wheel wheel, float theta, float h_s) 
{
    float theta_f = getKineticContactExitAngle(s,wheel,h_s);
    float theta_r = getKineticContactEntryAngle(s, wheel, h_s);
    float theta_m = (s.a_0 + (s.a_1 * wheel.s)) * theta_f;
    float k = pow(wheel.r, s.n) * ((s.k_c / wheel.b) + s.k_phi);
    float sigma = -1;

    if(theta_m < theta && theta < theta_f) 
    {
        sigma = k * pow(((cos(theta) - cos(theta_f))), s.n);
    } else if(theta_r < theta && theta < theta_m ) 
    {
        sigma = k * pow((cos(theta_f - (((theta - theta_r) / (theta_m - theta_r)) * (theta_f - theta_m))) - cos(theta_f)),s.n);
    }
    return sigma;
}

double BasicForceSolver::getStaticContactIntegrand(SoilPatch s, float theta, float theta_s) {
    return cos(theta) - (pow(cos(theta_s), s.n) * cos(theta)); 
}

float BasicForceSolver::getKineticContactExitAngle(SoilPatch s, Wheel wheel, float h_s) 
{
    return acos(1 - (h_s / wheel.r));   
}

float BasicForceSolver::getKineticContactEntryAngle(SoilPatch s, Wheel wheel, float h_s) 
{
    return -acos(1 - (s.l_s * h_s / wheel.r));
}

float BasicForceSolver::getTAUX(SoilPatch s, Wheel wheel, float theta, float h_s) 
{
    return (s.c + (getSigma(s,wheel,theta,h_s) * tan(s.phi))) * (1 - exp(-getJX(s,wheel,theta, h_s) / s.k_x));
}

float BasicForceSolver::getJX(SoilPatch s, Wheel wheel, float theta, float h_s) 
{
    return wheel.r * (getKineticContactExitAngle(s,wheel,h_s) - theta - ((1 - wheel.s) * (sin(getKineticContactExitAngle(s,wheel,h_s)) - sin(theta))));
}