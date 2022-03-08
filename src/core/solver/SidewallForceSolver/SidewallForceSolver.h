#ifndef _BASIC_FORCE_SOLVER
#define _BASIC_FORCE_SOLVER

#include "../WheelForceSolver.h"

namespace Hina {
    class SidewallForceSolver : public WheelForceSolver {

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
        SidewallForceSolver() {};
        ~SidewallForceSolver() {};

        void step();
        float getDrawbarPull();
        float getWheelSinkage();

        float getSidewallResistance(float beta, float theta_f) 


        float getDynamicSinkage(SoilPatch s, Wheel wheel);
        float getStaticSinkage(SoilPatch s, Wheel wheel);
        float getStaticContactAngle(SoilPatch s, Wheel wheel) const;
        float getKineticContactExitAngle(SoilPatch s, Wheel wheel);
        float getKineticContactEntryAngle(SoilPatch s, Wheel wheel);
        float getSigma(SoilPatch s, Wheel wheel, float theta);
        // float getJX(SoilPatch s, Wheel wheel, float theta);
        // float getTAUX(SoilPatch s, Wheel wheel, float theta);
        // float getFX(SoilPatch s, Wheel wheel, float theta, float theta_r, float theta_f);
        // float getFZ(SoilPatch s, Wheel wheel, float theta, float theta_r, float theta_f);
    };
}

#endif 