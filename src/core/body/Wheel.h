#ifndef WHEEL_H
#define WHEEL_H

#include "Body.h"

namespace Hina {
    class Wheel : Body {
    public:
        float W;
        float b;
        float beta;
        float r;
        float s;
    };
}

#endif