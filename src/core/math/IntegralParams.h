#ifndef INTEGRAL_PARAMS_H
#define INTEGRAL_PARAMS_H

#include "Function.h"
#include "Integral.h"

class Integral;

struct IntegralParams {
    Integral* integral = nullptr;
    Function integrand = nullptr;
    void* params = nullptr;
    double a;
    double b;
};

#endif