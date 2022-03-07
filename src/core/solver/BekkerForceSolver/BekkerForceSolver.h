#ifndef _BASIC_FORCE_SOLVER
#define _BASIC_FORCE_SOLVER

#include "../WheelForceSolver.h"

namespace Hina {
    class BekkerForceSolver : public WheelForceSolver {

    private:
        struct static_contact_integrand_param {
            SoilPatch s;
            double theta_s;
        };
        struct static_contact_integral_param {
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

        float getDynamicSinkage(SoilPatch s, Wheel wheel);
        float getStaticSinkage(SoilPatch s, Wheel wheel);
        float getStaticContactAngle(SoilPatch s, Wheel wheel) const;
        float getKineticContactExitAngle(SoilPatch s, Wheel wheel, float h_k);
        float getKineticContactEntryAngle(SoilPatch s, Wheel wheel, float h_k);
        float getSigma(SoilPatch s, Wheel wheel, float theta, float h_k);
        float getJX(SoilPatch s, Wheel wheel, float theta, float h_k);
        float getTAUX(SoilPatch s, Wheel wheel, float theta, float h_k);
        float getFX(SoilPatch s, Wheel wheel, float theta, float theta_r, float theta_f);
        float getFZ(SoilPatch s, Wheel wheel, float theta, float theta_r, float theta_f);
    };
}

#endif 