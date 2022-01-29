#include "BekkerForceSolver.h"
#include <gsl/gsl_integration.h>
#include <functional>
#include <iostream>
#include <math.h>

using namespace Hina;
using namespace std;
using namespace std::placeholders; 

void BekkerForceSolver::step() 
{

}

float BekkerForceSolver::getDrawbarPull() 
{
    return 0;
}   

/*
 * Get the kinematic sinkage of the wheel.
 * 
 */
float BekkerForceSolver::getWheelSinkage() 
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
    return 0;
}

/*
 * Calculates the static contact angle of the wheel.
 *
 * @param s The soil patch.
 * @param wheel The wheel.
 * @return The static contact angle of the wheel.
 */
float BekkerForceSolver::getStaticContactAngle(SoilPatch s, Wheel wheel) 
{
    float k = (wheel.r, s.n + 1) * (s.k_c + (s.k_phi * wheel.b));
    static_integrand_param* p = new static_integrand_param();
    p->s = s;
    p->theta_s = 0.9;

    auto integrand = [](double theta, void * params) -> double {
        auto p = (struct static_integrand_param *) params;
        return cos(theta) - (pow(cos(p->theta_s), p->s.n) * cos(theta)); 
    };

    auto integral = [this, &k, &s, &integrand, &p] (float theta_s) -> double {
        p->theta_s = theta_s;
        
        double result, error;
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_function F;
        F.function = integrand;
        F.params = &p;
        gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);
        printf("result = %.18f\n", result);
        gsl_integration_workspace_free (w);

        return k * result;
    };
    std::cout << integral(0.4) << std::endl;

    delete (p);
    return -1;
}

float BekkerForceSolver::getSigma(SoilPatch s, Wheel wheel, float theta, float h_s) 
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

double BekkerForceSolver::getStaticContactIntegrand(SoilPatch s, float theta, float theta_s) {
    return cos(theta) - (pow(cos(theta_s), s.n) * cos(theta)); 
}

float BekkerForceSolver::getKineticContactExitAngle(SoilPatch s, Wheel wheel, float h_s) 
{
    return acos(1 - (h_s / wheel.r));   
}

float BekkerForceSolver::getKineticContactEntryAngle(SoilPatch s, Wheel wheel, float h_s) 
{
    return -acos(1 - (s.l_s * h_s / wheel.r));
}

float BekkerForceSolver::getTAUX(SoilPatch s, Wheel wheel, float theta, float h_s) 
{
    return (s.c + (getSigma(s,wheel,theta,h_s) * tan(s.phi))) * (1 - exp(-getJX(s,wheel,theta, h_s) / s.k_x));
}

float BekkerForceSolver::getJX(SoilPatch s, Wheel wheel, float theta, float h_s) 
{
    return wheel.r * (getKineticContactExitAngle(s,wheel,h_s) - theta - ((1 - wheel.s) * (sin(getKineticContactExitAngle(s,wheel,h_s)) - sin(theta))));
}