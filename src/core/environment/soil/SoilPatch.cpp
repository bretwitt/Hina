#include "SoilPatch.h"

using namespace Hina;

SoilPatch::SoilPatch() {

}

SoilPatch::~SoilPatch() {

}

void SoilPatch::init_parameters(float c, float phi, float k_c, float k_phi, float n, float a_0, float a_1, float rho_d, float l_s, float k_x, float k_y) {
    this->c = c;
    this->phi = phi;
    this->k_c = k_c;
    this->k_phi = k_phi;
    this->n = n;
    this->a_0 = a_0;
    this->a_1 = a_1;
    this->rho_d = rho_d;
    this->l_s = l_s;
    this->k_x = k_x;
    this->k_y = k_y;
}