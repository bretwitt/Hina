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
 * 
 */
float BekkerForceSolver::getDynamicSinkage(SoilPatch s, Wheel w) 
{
    struct integrand_params {
        BekkerForceSolver* solver;
        double h;
        SoilPatch s;
        Wheel w;
    };
    struct integral_params {
        BekkerForceSolver* solver;
        SoilPatch s;
        Wheel w;
        gsl_func_ptr integrand;
    };
    struct function_params {
        BekkerForceSolver* solver;
        Wheel wheel;
        SoilPatch s;
        gsl_func_ptr integral;
        gsl_func_ptr integrand;
    };

    auto integrand = [](double theta, void* params) -> double {
        auto p = *(struct integrand_params*) params;
        return (p.solver->getTAUX(p.s,p.w,theta) * sin(theta)) + (p.solver->getSigma(p.s,p.w,theta) * cos(theta));
    };

    auto integral = [](double h, void* params) -> double {
        auto p = *(struct integral_params*) params;

        auto p1 = new integrand_params {
            .solver = p.solver,
            .h = h,
            .s = p.s,
            .w = p.w
        };

        double theta_r = p.solver->getKineticContactEntryAngle(p.s,p.w);
        double theta_f = p.solver->getKineticContactExitAngle(p.s,p.w);

        double result, error;
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

        gsl_function F;
        F.function = p.integrand;
        F.params = p1;

        gsl_integration_qags (&F, theta_r, theta_f, 0, 1e-7, 1000, w, &result, &error);

        gsl_integration_workspace_free (w);
        delete(p1);

        return result;
    };

    auto fn = [](double h, void* params) -> double {
        auto p = *(struct function_params *) params;
        auto integr_params = new integral_params {
            .solver = p.solver,            
            .s = p.s,
            .w = p.wheel,
            .integrand = p.integrand
        };
        double res = p.integral(h, integr_params);
        delete(integr_params);
        return p.wheel.W - (p.wheel.r * p.wheel.b * res);
    };

    auto p = new function_params {
        .solver = this,
        .wheel = wheel,
        .s = s,
        .integral = integral,
        .integrand = integrand
    };

    auto params = new integral_params {
        .solver = this,
        .s = s,
        .w = w,
        .integrand = integrand,
    };

    // auto params = new integrand_params {
    //     .solver = this,
    //     .h = 0.5,
    //     .s = s,
    //     .w = w
    // };

    vector<double> x,y;

    for(double i = -0.783f; i < 0.645f; i += 0.01) {
        double res = integral(i, params);
        y.push_back(res);
        x.push_back(i);
    }

    plt::plot(x,y);
    plt::show();

    delete(params);
    // Solve for h s.t fn(h) = 0 
    // double res = fn(0.5f, (void*)(p));
    double res = -1;
    
    delete p;

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
        gsl_integration_qags (&F, -theta_s, theta_s, 0, 1e-7, 1000, w, &result, &error);
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

float BekkerForceSolver::getSigma(SoilPatch s, Wheel wheel, float theta) 
{
    float h_s = getStaticSinkage(s,wheel);
    float theta_f = getKineticContactExitAngle(s,wheel);
    float theta_r = getKineticContactEntryAngle(s, wheel);
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

float BekkerForceSolver::getTAUX(SoilPatch s, Wheel wheel, float theta) 
{
    // std::cout <<  (1 - exp(-getJX(s,wheel,theta) / s.k_x)) << std::endl;
    return (s.c + (getSigma(s,wheel,theta) * tan(s.phi))) * (1 - exp(-getJX(s,wheel,theta) / s.k_x));
}

float BekkerForceSolver::getJX(SoilPatch s, Wheel wheel, float theta) 
{
    float h_s = getStaticSinkage(s,wheel);
    return wheel.r * (getKineticContactExitAngle(s,wheel) - theta - ((1 - wheel.s) * (sin(getKineticContactExitAngle(s,wheel)) - sin(theta))));
}

