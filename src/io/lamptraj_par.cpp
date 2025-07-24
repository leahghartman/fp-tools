#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include "core/base.h"
#include "io/lamptraj_par.h"

LAMMPSTrajectory::LAMMPSTrajectory(const std::string& filename) {
  file.open(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open LAMMPS file: " + filename);
  }
}

/*
* Reads the next frame from a LAMMPS trajectory dump file and populates the given Frame object.
* 
* This function scans through the file to locate the start of the next frame by looking for the
* "ITEM: TIMESTEP" marker. Once found, it reads and parses:
*   - The simulation timestep
*   - Number of atoms
*   - Box boundaries for each spatial dimension
*   - Atom column
*   - Atom data for each particle in the frame
*
*/
bool LAMMPSTrajectory::next_frame(Frame& frame) {
  std::string line;

  // Seek frame start
  while (std::getline(file, line)) {
    if (line.find("ITEM: TIMESTEP") != std::string::npos) break;
  }
  if (file.eof()) return false;

  // TIMESTEP
  std::getline(file, line);
  frame.timestep = std::strtol(line.c_str(), nullptr, 10);

  // NUMBER OF ATOMS
  std::getline(file, line);  // skip header
  std::getline(file, line);
  frame.atomnum = std::strtol(line.c_str(), nullptr, 10);

  // BOX BOUNDS
  std::getline(file, line);  // skip header
  for (int i = 0; i < 3; ++i) {
    std::getline(file, line);
    const char* ptr = line.c_str();
    char* end;
    frame.box_bounds[i][0] = std::strtod(ptr, &end);
    frame.box_bounds[i][1] = std::strtod(end, nullptr);
  }

  // ATOM HEADER
  std::getline(file, line);
  std::unordered_map<std::string, int> column_index;  // Map: column name → index

  {
    std::istringstream iss(line.substr(12));  // skip "ITEM: ATOMS "
    std::string col;
    int idx = 0;
    while (iss >> col) {
      column_index[col] = idx++;
    }
  }

  // Ensure required columns exist
  bool has_id = column_index.count("id") > 0;
  if (!has_id) {
    std::cout << "Warning: 'id' column missing — assigning sequential atom IDs.\n";
  }
  int id_idx = has_id ? column_index["id"] : -1;
  if (!(
        column_index.count("xu") || column_index.count("xsu") ||
        column_index.count("xs") || column_index.count("x")
      )) throw std::runtime_error("Missing x coordinate column");

  if (!(
        column_index.count("yu") || column_index.count("ysu") ||
        column_index.count("ys") || column_index.count("y")
      )) throw std::runtime_error("Missing y coordinate column");

  if (!(
        column_index.count("zu") || column_index.count("zsu") ||
        column_index.count("zs") || column_index.count("z")
      )) throw std::runtime_error("Missing z coordinate column");

  // Resolve coordinate column names
  auto find_coord_index = [&](const std::string& base) -> int {
    if (column_index.count(base + "u"))  return column_index[base + "u"];
    if (column_index.count(base + "su")) return column_index[base + "su"];
    if (column_index.count(base + "s"))  return column_index[base + "s"];
    if (column_index.count(base))        return column_index[base];
    throw std::runtime_error("Missing coordinate column for: " + base);
  };

  int type_idx = -1;
  if (column_index.count("type")) {
    type_idx = column_index["type"];
  } else if (column_index.count("element")) {
    type_idx = column_index["element"];
  } else {
    // Could add more fallbacks or keep -1
    std::cout << "Warning: no 'type' or 'element' column found, atom types will be 'unknown'.\n";
  }
  int x_idx    = find_coord_index("x");
  int y_idx    = find_coord_index("y");
  int z_idx    = find_coord_index("z");

  int vx_idx = column_index.count("vx") ? column_index["vx"] : -1;
  int vy_idx = column_index.count("vy") ? column_index["vy"] : -1;
  int vz_idx = column_index.count("vz") ? column_index["vz"] : -1;

  bool has_velocity = (vx_idx >= 0 && vy_idx >= 0 && vz_idx >= 0);


  // ATOM DATA
  frame.atoms.clear();
  frame.atoms.reserve(frame.atomnum);

  std::vector<const char*> token_ptrs;
  std::vector<size_t> token_lens;
  token_ptrs.reserve(16);
  token_lens.reserve(16);

  for (int i = 0; i < frame.atomnum; ++i) {
    std::getline(file, line);
    token_ptrs.clear();
    token_lens.clear();

    const char* ptr = line.c_str();
    while (*ptr) {
      while (*ptr == ' ' || *ptr == '\t') ++ptr;
      if (!*ptr) break;
      const char* start = ptr;
      while (*ptr != ' ' && *ptr != '\t' && *ptr != '\0') ++ptr;
      token_ptrs.push_back(start);
      token_lens.push_back(ptr - start);
    }

    Atom atom;
    atom.id = has_id ? std::strtol(token_ptrs[id_idx], nullptr, 10) : (i + 1);
    if (type_idx >= 0) {
      atom.type = std::string(token_ptrs[type_idx], token_lens[type_idx]);
    } else {
      atom.type = "unknown";  // or "" or std::nullopt if optional
    }
    atom.x = std::strtod(token_ptrs[x_idx], nullptr);
    atom.y = std::strtod(token_ptrs[y_idx], nullptr);
    atom.z = std::strtod(token_ptrs[z_idx], nullptr);

    if (has_velocity) {
      atom.vx = std::strtod(token_ptrs[vx_idx], nullptr);
      atom.vy = std::strtod(token_ptrs[vy_idx], nullptr);
      atom.vz = std::strtod(token_ptrs[vz_idx], nullptr);
    } else {
      atom.vx = atom.vy = atom.vz = 0.0;
    }
    frame.atoms.push_back(atom);
  }

  return true;
}
