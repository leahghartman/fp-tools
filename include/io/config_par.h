/*
* -----------------------------------------------------------------------------
* Configuration Structures
* -------------------------
* Defines the data structures used to represent the user input configuration
* for trajectory analysis. Supports input, output, and various analysis modules
* including MSD, RDF, FKT, FSKT, MFPT, and D(r).
*
* These structs are designed to be populated by parsing a YAML configuration file.
*
* Authors: leahghartman
* -----------------------------------------------------------------------------
*/

#pragma once
#include <string>
#include <vector>

// --- Input Configuration ---
struct InputConfig {
    std::string file{};           // Path to the input trajectory file; required
    std::string coord_type{};     // Coordinate type, e.g. "wrapped" or "unwrapped"; optional
    int start_frame{0};           // Starting frame index (inclusive); default 0
    int end_frame{-1};            // Ending frame index (inclusive); default -1 (last frame) 
    int frame_interval{1};        // Interval to sample frames, e.g. 1 means every frame
};

// --- Output Configuration ---
struct OutputConfig {
  std::string path{"output"}; 
  int verbosity{1};          // Verbosity level: 0 = silent, 1 = normal
  int status_interval{10};   // Interval (in frames) to print status/progress updates
};

// --- MSD Configuration ---
struct MSDConfig {
    int max_lag{-1};            // Maximum lag time (frames) for MSD calculation; -1 means full length
    std::string plot_format{};  // Format for output plot files, e.g. "png"; empty disables plotting
};

// --- Container for ALL analysis configurations ---
struct AnalysisConfig {
    MSDConfig msd{};    // MSD parameters
};

// --- Top-level configurationi struct aggregating all sections ---
struct Config {
    InputConfig input{};        // Input file and frame selection options
    OutputConfig output{};      // Output verbosity and status options
    AnalysisConfig analysis{};  // Analysis-specific parameters
};

// Parser function that reads from a YAML file and returns a populated Config struct
Config parse_config(const std::string& filename);
