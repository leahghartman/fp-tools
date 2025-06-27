#pragma once
#include <vector>           // Used in LammpsFrame.atoms
#include <string>           // For filename input
#include <unordered_map>    // For column indexing during parsing

// Define a data structure for each atom
struct Atom {
  int id{};                   // Atom ID (unique identifier in LAMMPS)
  std::string type{};                 // Atom type (integer or atomic symbol, e.g., "1", "Al")
  double x{}, y{}, z{};       // Coordinates in 3D space
};

// Define a data structure for each trajectory frame (variable names are based on the format of the dump files)
struct Frame {
  int timestep{};                     // Simulation timestep number
  int atomnum{};
  double box_bounds[3][2]; 
  std::vector<Atom> atoms{};          // An array of Atom objects in a frame; initialized to be empty
  std::string coord_type{};           // "unscaled", "scaled", "unwrapped", "scaled_unwrapped"
};

struct FrameCollection {
  std::vector<Frame> frames;
  std::string coord_type;  // "unscaled", "scaled", "unwrapped", or "scaled_unwrapped"
};
