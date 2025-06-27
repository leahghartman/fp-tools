//  TODO: Finish documentation

/* 
* -----------------------------------------------------------------------------
* XYZ Trajectory Parser
* ----------------------
* This file implements the XYZTrajectory class for reading VASP XDACAR-style
* trajectory files. It supports only rectangular (orthorhombic) boxes and
* assumes atomic positions are provided in fractional ("Direct") format.
*
* The parser reads the box geometry, element types, atom counts, and then
* returns one Frame at a time, scaled into real coordinates.
*
* Authors: leahghartman
* -----------------------------------------------------------------------------
*/

#include <sstream>
#include "io/xyztraj_par.h"

XYZTrajectory::XYZTrajectory(const std::string& filename) {
    file.open(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open XYZ file: " + filename);
    }
}

bool XYZTrajectory::next_frame(Frame& frame) {
    std::string line{};
    int natoms{};

    // Read number of atoms
    if (!std::getline(file, line)) return false;  // End of file
    std::istringstream(line) >> natoms;

    // Skip comment line
    if (!std::getline(file, line)) return false;

    frame.atoms.clear();
    frame.atoms.reserve(natoms);
    frame.timestep = timestep_counter++;
    frame.atomnum = natoms;

    // Parse atoms
    for (int i{0}; i < natoms; ++i) {
        if (!std::getline(file, line)) return false;

        std::istringstream iss(line);
        std::string element;
        double x, y, z;

        if (!(iss >> element >> x >> y >> z)) {
            throw std::runtime_error("Invalid atom line in XYZ file");
        }

        Atom atom;
        atom.type = "";  // Optional: convert `element` to int if needed
        atom.x = x;
        atom.y = y;
        atom.z = z;

        frame.atoms.push_back(atom);
    }
    return true;
}
