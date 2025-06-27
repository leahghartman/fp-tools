#pragma once
#include "core/base.h"

class TrajectoryProvider {
public:
    virtual ~TrajectoryProvider() = default;

    // Pure virtual method to be implemented by derived classes
    virtual bool next_frame(Frame& frame) = 0;
};

