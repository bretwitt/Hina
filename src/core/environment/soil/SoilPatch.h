#ifndef EE061387_B93D_4CA5_8939_8CFD3A8A955A
#define EE061387_B93D_4CA5_8939_8CFD3A8A955A

namespace Hina {
    class SoilPatch {
        
    public:
        SoilPatch();
        ~SoilPatch();
        void init_parameters(float c, float phi, float k_c, float k_phi, float n, float a_0, float a_1, float rho_d, float l_s, float k_x, float k_y);
        float c;
        float phi;
        float k_c;
        float k_phi;
        float n;
        float a_0;
        float a_1;
        float rho_d;
        float l_s;
        float k_x;
        float k_y;
    };
}

#endif /* EE061387_B93D_4CA5_8939_8CFD3A8A955A */
