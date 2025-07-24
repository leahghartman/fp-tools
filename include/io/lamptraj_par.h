#pragma once
#include <string>
#include <fstream>
#include "core/base.h"
#include "core/traj_prov.h"

class LAMMPSTrajectory : public TrajectoryProvider {
public:
  explicit LAMMPSTrajectory(const std::string& filename);
  bool next_frame(Frame& frame) override;

private:
  std::ifstream file;
  std::string coord_type{};
};
