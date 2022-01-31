#ifndef _BASIC_FORCE_SOLVER
#define _BASIC_FORCE_SOLVER

#include "../WheelForceSolver.h"

namespace Hina {
    class BekkerForceSolver : public WheelForceSolver {

    private:
        using integrand_ptr = double(*)(double, void*);
        struct static_integrand_param {
            SoilPatch s;
            double theta_s;
        };
        struct integral_param {
            double k;
            SoilPatch s;
            Wheel wheel;
            double (*integrand) (double, void *);
        };

    public:
        BekkerForceSolver() {};
        ~BekkerForceSolver() {};

        void step();
        float getDrawbarPull();
        float getWheelSinkage();

        float getStaticSinkage(SoilPatch s, Wheel wheel);
        float getStaticContactAngle(SoilPatch s, Wheel wheel);
        double getStaticContactIntegrand(SoilPatch s, float theta, float theta_s);
        float getKineticContactExitAngle(SoilPatch s, Wheel wheel, float h_s);
        float getKineticContactEntryAngle(SoilPatch s, Wheel wheel, float h_s);
        float getSigma(SoilPatch s, Wheel wheel, float theta, float h_s);
        float getJX(SoilPatch s, Wheel wheel, float theta, float h_s);
        float getTAUX(SoilPatch s, Wheel wheel, float theta, float h_s);
        float getFX(SoilPatch s, Wheel wheel, float theta, float theta_r, float theta_f);
        float getFZ(SoilPatch s, Wheel wheel, float theta, float theta_r, float theta_f);
    };
}

#endif 