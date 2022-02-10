#ifndef _INTEGRAL_H
#define _INTEGRAL_H

#include "IntegralParams.h"
#include "Function.h"
#include <gsl/gsl_integration.h>

class Integral : public Function {


public:
    function_ptr getFunction() {
        return [] (double n, void* params) -> double {
            auto p = *(struct IntegralParams*) params;
            return p.integral->evaluate(p);
        };
    }

    double evaluate(IntegralParams params) {
        double result, error;
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

        gsl_function F;
        F.function = params.integrand.getFunction();
        F.params = params.params;

        gsl_integration_qags (&F, params.a, params.b, 0, 1e-7, 10, w, &result, &error);
        gsl_integration_workspace_free (w);
        return result;
    };

};

#endif
