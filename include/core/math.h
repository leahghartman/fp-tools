#pragma once

#include <array>
#include <vector>
#include <cstddef>
#include "core/base.h"

// Vector Operations
inline double squared_displacement(const Atom& a, const Atom& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}
