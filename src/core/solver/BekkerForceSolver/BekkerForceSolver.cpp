#include "BekkerForceSolver.h"
#include <gsl/gsl_integration.h>
#include <functional>
#include <iostream>
#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace Hina;
using namespace std;

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
    float res = -1;

    auto p = new static_integrand_param();
    p->s = s;
    p->theta_s = 0;

    auto integrand = [](double theta, void * params) -> double {
        auto p = (struct static_integrand_param *) params;
        return cos(theta) - (pow(cos(p->theta_s), p->s.n) * cos(theta)); 
    };

    auto integral = [] (double theta_s, void* params) -> double {
        auto p = (struct integral_param *) params;
        auto p1 = new static_integrand_param {
            .s = p->s,
            .theta_s = theta_s
        };

        double result, error;
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_function F;
        F.function = p->integrand;
        F.params = p1;
        gsl_integration_qags (&F, -theta_s, theta_s, 0, 1e-7, 1000, w, &result, &error);
        gsl_integration_workspace_free (w);
        return p->wheel.W - (p->k * result);
    };

    auto integral_param_struct = new integral_param {
        .k = k,
        .p = p,
        .s = s,
        .wheel = wheel,
        .integrand = integrand
    };

    std::cout << "Starting static contact angle search..." << std::endl;

    gsl_function F;
    F.function = integral;
    F.params = integral_param_struct;
    gsl_root_fsolver * sol = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
    gsl_root_fsolver_set (sol, &F, 0, M_PI_2);
    int status;
    int iter = 0;
    do {
        iter++;
        status = gsl_root_fsolver_iterate (sol);
        res = gsl_root_fsolver_root (sol);
        status = gsl_root_test_residual (res, 1e-7);
        std::cout << "iter: " << iter << " res: " << res << " status: " << status << std::endl;
    } while (status == GSL_CONTINUE && iter < 1000);
    
    gsl_root_fsolver_free (sol);

    delete (p);
    delete (integral_param_struct);

    return res;
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