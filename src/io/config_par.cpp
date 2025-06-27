//  TODO: Finish documentation
//  TODO: Finish adding inputs that Jerome/Alfred want

/* 
* -----------------------------------------------------------------------------
* Configuration Parser
* ---------------------
* This file implements the VASPTrajectory class for reading VASP XDATCAR-style
* trajectory files. It supports only rectangular (orthorhombic) boxes and
* assumes atomic positions are provided in fractional ("Direct") format.
*
* The parser reads the box geometry, element types, atom counts, and then
* returns one Frame at a time, scaled into real coordinates.
*
* Authors: leahghartman
* -----------------------------------------------------------------------------
*/

#include <stdexcept>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include "io/config_par.h"

Config parse_config(const std::string& filename) {
  Config config{};

  // Parse the YAML file
  YAML::Node root = YAML::LoadFile(filename);

  // --- [input] section (required) ---
  const auto& input{root["input"]};
  config.input.file = input["file"].as<std::string>();

  config.input.start_frame = input["start_frame"] ? input["start_frame"].as<int>() : 0;
  config.input.end_frame = input["end_frame"] ? input["end_frame"].as<int>() : -1;
  config.input.frame_interval = input["frame_interval"] ? input["frame_interval"].as<int>() : 1;

  // --- [output] section (optional) ---
  if (root["output"]) {
    const auto& output{root["output"]};
    config.output.path = output["path"] ? output["path"].as<std::string>() : "output";
    config.output.verbosity = output["verbosity"] ? output["verbosity"].as<int>() : 1;
    config.output.status_interval = output["status_interval"] ? output["status_interval"].as<int>() : 10;
  }

  // --- [analysis] (optional) ---
  if (root["analysis"]) {
    const auto& analysis = root["analysis"];

    // MSD Configuration
    if (analysis["msd"]) {
      const auto& msd{analysis["msd"]};
      config.analysis.msd.max_lag = msd["max_lag"] ? msd["max_lag"].as<int>() : -1;
      config.analysis.msd.plot_format = msd["plot_format"] ? msd["plot_format"].as<std::string>() : "";
    }
  }
  return config;
}
