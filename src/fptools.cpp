/*
 * -----------------------------------------------------------------------------
 * fptools.cpp
 * -----------------------------------------------------------------------------
 * This is the main entry point for fp-tools, a C++ molecular dynamics analysis
 * library supporting analysis techniques like MSD, RDF, FKT, FSKT, MFPT, and D(r).
 *
 * This file:
 *   - Loads and parses the YAML configuration file
 *   - Detects the trajectory format (e.g. LAMMPS, VASP, XYZ)
 *   - Iterates over frames from the trajectory file
 *   - Calls selected analysis routines based on user config
 *   - Optionally creates plots via external Python scripts
 *
 * Authors: leahghartman
 * -----------------------------------------------------------------------------
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <string_view>
#include <yaml-cpp/yaml.h>

// Core data structures and config parsing
#include "core/traj_prov.h"
#include "io/config_par.h"

// Input trajectory parsers
//#include "io/lamptraj_par.h"
//#include "io/vasptraj_par.h"
#include "io/xyztraj_par.h"

// Analysis modules
#include "analysis/msd.h"
//#include "analysis/rdf.h"

namespace fs = std::filesystem;
fs::path get_executable_dir(char* argv0) {
    // This returns the directory where the executable resides as an absolute path
    return fs::weakly_canonical(fs::absolute(argv0)).parent_path();
}

/*
 * -----------------------------------------------------------------------------
 * Helper function to determine trajectory format from filename.
 *
 * This infers file type as there is currently non input for the user to 
 * explicitly specify it. For now, returns lowercase version of extension 
 * or fallback to a blank string, which throws an error in main().
 * -----------------------------------------------------------------------------
 */
std::string infer_format(const std::string& filename) {
    //if (filename.find("XDATCAR") != std::string::npos) return "vasp";
    if (filename.ends_with(".xyz")) return "xyz";
    //if (filename.ends_with(".lammpstrj")) return "lammps";
    return "";  // Default fallback
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: fptools_main <config.yaml>" << std::endl;
    return 1;
  }

  std::string config_file{argv[1]};

  // ---------------------------------------------
  // Load and parse configuration file
  // ---------------------------------------------
  Config config = parse_config(config_file);
  std::string format{infer_format(config.input.file)};
  std::string output_dir{config.output.path};

  // ---------------------------------------------
  // Parse the trajectory/dump file
  // ---------------------------------------------

  // Load trajectory based on file extension or config.input.coord_type (simplified here)
  std::unique_ptr<TrajectoryProvider> traj{};

  // Simple heuristic by extension for example; ideally use config.input.format or coord_type
  //if (config.input.file.find(".lammpstrj") != std::string::npos) {
    //traj = std::make_unique<LAMMPSTrajectory>{config.input.file};
  //} else if (config.input.file.find("XDATCAR") != std::string::npos) {
    //traj = std::make_unique<VASPTrajectory>{config.input.file};
  if (config.input.file.find(".xyz") != std::string::npos) {
    traj = std::make_unique<XYZTrajectory>(config.input.file);
  } else {
    std::cerr << "Unsupported trajectory format or file extension.\n";
    return 1;
  }

  // Collect frames according to start, end, interval
  std::vector<Frame> frames{};
  int frame_count{0};
  Frame frame{};
  while (traj->next_frame(frame)) {
    if (frame_count >= config.input.start_frame &&
      (config.input.end_frame < 0 || frame_count <= config.input.end_frame) &&
      ((frame_count - config.input.start_frame) % config.input.frame_interval == 0)) {
      frames.push_back(frame);
    }
    ++frame_count;
  }

  if (frames.empty()) {
    std::cerr << "No frames loaded based on start/end/frame_interval settings." << std::endl;
    return 1;
  }

  // Verbosity and status interval for output
  int verbosity{config.output.verbosity};
  int status_interval{config.output.status_interval};

  // ---------------------------------------------
  // Perform selected analysis
  // ---------------------------------------------
  
  try {
    std::cout << "Current working directory: " << fs::current_path() << "\n";
    std::cout << "Requested output directory: " << output_dir << "\n";

    if (fs::exists(output_dir)) {
        std::cout << "Output directory already exists.\n";
    } else {
        bool created = fs::create_directories(output_dir);
        std::cout << "Directory creation attempt: " << (created ? "Success" : "Failed") << "\n";

        if (!created) {
            std::cerr << "Error: Could not create output directory: " << output_dir << "\n";
            return 1;
        }
    }
  } catch (const std::exception& e) {
    std::cerr << "Exception caught while creating directory: " << e.what() << "\n";
    return 1;
  }

  // ---- MSD ----
  if (config.analysis.msd.max_lag > 0) {
    if (verbosity > 0) std::cout << "Computing MSD...\n";

    MSDAccumulator msd_accum(config.analysis.msd.max_lag);
    std::string output_file{output_dir + "msd.dat"};
    
    msd_accum = compute_msd(frames, config.input.start_frame, config.input.end_frame, config.analysis.msd.max_lag);
    msd_accum.write(output_file);

    if (!config.analysis.msd.plot_format.empty()) {
      msd_accum.plot(output_file, config.analysis.msd.plot_format);
    }
  }
  if (verbosity > 0) std::cout << "All requested analyses completed.\n";
  return 0;
}
