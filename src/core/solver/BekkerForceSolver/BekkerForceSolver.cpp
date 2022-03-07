#include "BekkerForceSolver.h"
#include <gsl/gsl_integration.h>
#include <functional>
#include <iostream>
#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;
using namespace Hina;
using namespace std;

void BekkerForceSolver::step() 
{

}

float BekkerForceSolver::getDrawbarPull() 
{
    return 0;
}   

float BekkerForceSolver::getWheelSinkage() {
    return getStaticSinkage(this->soil,this->wheel) + getDynamicSinkage(this->soil,this->wheel);
}

/*
 * Get the kinematic sinkage of the wheel.
 */
float BekkerForceSolver::getDynamicSinkage(SoilPatch s, Wheel w) 
{
    struct integrand_params {
        BekkerForceSolver* solver;
        double h_k;
        SoilPatch s;
        Wheel w;
    };

    struct integral_params {
        double (*integrand) (double, void *);
        BekkerForceSolver* solver;
        SoilPatch s;
        Wheel w;
    };

    auto integrand = [](double theta, void* params) -> double {
        auto p = *(struct integrand_params*) params;
        double res = (p.solver->getTAUX(p.s, p.w, theta, p.h_k) * sin(theta)) + (p.solver->getSigma(p.s, p.w, theta, p.h_k) * cos(theta));
        return res;
    };

    auto integral = [](double h_k, void* params) -> double {
        auto p = *(struct integral_params*) params;
        auto i_p = new integrand_params { 
            .solver = p.solver,
            .h_k = h_k,
            .s = p.s,
            .w = p.w
        };
        std::cout << h_k << std::endl;

        double quad;
        double error;

        double theta_f = p.solver->getKineticContactEntryAngle(p.s, p.w, h_k);
        double theta_r = p.solver->getKineticContactExitAngle(p.s, p.w, h_k);

        std::cout << "theta_f " << theta_f << std::endl;
        std::cout << "theta_r " << theta_r << std::endl;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_function F;
        F.function = p.integrand;
        F.params = i_p;
        gsl_integration_qags (&F, theta_r, theta_f, 0, 1e-7, 200, w, &quad, &error);
        gsl_integration_workspace_free (w);

        delete i_p;

        return p.w.W - (p.w.r * p.w.b * quad);
    };

    auto params = new integral_params {
        .integrand = integrand,
        .solver = this,
        .s = s,
        .w = w
    };

    double res;
    gsl_function F;
    F.function = integral;
    F.params = params;
    gsl_root_fsolver * sol = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
    gsl_root_fsolver_set (sol, &F, 0.001, 0.2);
    int status;
    int iter = 0;
    do {
        iter++;
        status = gsl_root_fsolver_iterate (sol);
        res = gsl_root_fsolver_root (sol);
        status = gsl_root_test_residual (res, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);
    
    gsl_root_fsolver_free (sol);

    // std::cout << integral(0.016, params) << std::endl;
    return res;
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

/*
 * Calculates the static contact angle of the wheel.
 *
 * @param s The soil patch.
 * @param wheel The wheel.
 * @return The static contact angle of the wheel.
 */
float BekkerForceSolver::getStaticContactAngle(SoilPatch s, Wheel wheel) const
{
    float k = pow(wheel.r, s.n + 1) * (s.k_c + (s.k_phi * wheel.b)); 
      
    float res = -1;

    auto integrand = [](double theta, void * params) -> double {
        auto p = (struct static_contact_integrand_param *) params;
        return pow(cos(theta) - cos(p->theta_s), p->s.n) * cos(theta); 
    };

    auto integral = [] (double theta_s, void* params) -> double {
        auto p = (struct static_contact_integral_param *) params;
        auto p1 = new static_contact_integrand_param {
            .s = p->s,
            .theta_s = theta_s
        };

        double result, error;
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_function F;
        F.function = p->integrand;
        F.params = p1;
        gsl_integration_qags (&F, -theta_s, theta_s, 0, 1e-7, 200, w, &result, &error);
        gsl_integration_workspace_free (w);

        double W = p->wheel.W;
        double k = p->k;

        delete(p1);
        return W - (k * result);
    };

    auto integral_param_struct = new static_contact_integral_param {
        .k = k,
        .s = s,
        .wheel = wheel,
        .integrand = integrand
    };

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
    } while (status == GSL_CONTINUE && iter < 1000);
    
    gsl_root_fsolver_free (sol);

    delete (integral_param_struct);

    return res;
}

float BekkerForceSolver::getSigma(SoilPatch s, Wheel wheel, float theta, float h_k) 
{
    float h_s = getStaticSinkage(s,wheel);
    float theta_f = getKineticContactEntryAngle(s, wheel, h_k);
    float theta_r = getKineticContactExitAngle(s, wheel, h_k);
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

float BekkerForceSolver::getKineticContactEntryAngle(SoilPatch s, Wheel wheel, float h_k) 
{
    return acos(1 - (h_k / wheel.r));   
}

float BekkerForceSolver::getKineticContactExitAngle(SoilPatch s, Wheel wheel, float h_k) 
{    
    return -acos(1 - (s.l_s * h_k / wheel.r));
}

float BekkerForceSolver::getTAUX(SoilPatch s, Wheel wheel, float theta, float h_k) 
{
    return (s.c + (getSigma(s,wheel,theta,h_k) * tan(s.phi))) * (1 - exp(-getJX(s, wheel, theta, h_k / s.k_x)));
}

float BekkerForceSolver::getJX(SoilPatch s, Wheel wheel, float theta, float h_k) 
{
    return wheel.r * (getKineticContactEntryAngle(s, wheel, h_k) - theta - ((1 - wheel.s) * (sin(getKineticContactEntryAngle(s, wheel, h_k)) - sin(theta))));
}

