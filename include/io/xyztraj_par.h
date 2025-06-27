#pragma once
#include <fstream>
#include <string>
#include "core/base.h"
#include "core/traj_prov.h"

class XYZTrajectory : public TrajectoryProvider {
public:
    explicit XYZTrajectory(const std::string& filename);

    // Override the pure virtual method from base class
    bool next_frame(Frame& frame) override;

private:
    std::ifstream file{};
    int timestep_counter{0};
};

