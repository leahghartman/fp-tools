#pragma once
#include <vector>           // Used in LammpsFrame.atoms
#include <string>           // For filename input
#include <unordered_map>    // For column indexing during parsing

struct Point3 {
    double x, y, z;
    int original_idx; // Storing original atom index is crucial
    int type_id;      // Storing type_id can speed up filtering
};

struct PointCloudAdapter {
    const std::vector<Point3>& pts;

    PointCloudAdapter(const std::vector<Point3>& p) : pts(p) {}

    inline size_t kdtree_get_dims() const { return 3; }
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // *** THIS IS THE CRUCIAL CHANGE FOR Error 3 ***
    // Removed template <class T> to be explicit about `double`
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return pts[idx].x;
        else if (dim == 1) return pts[idx].y;
        else return pts[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

// Define a data structure for each atom
struct Atom {
  int id{};                   // Atom ID (unique identifier in LAMMPS)
  std::string type{};         // Atom type (integer or atomic symbol, e.g., "1", "Al")
  std::string species{};
  double x{}, y{}, z{};       // Coordinates in 3D space
  double vx{}, vy{}, vz{};
  int molecule_id = -1;
};

// Define a data structure for each trajectory frame (variable names are based on the format of the dump files)
struct Frame {
  int timestep{};                     // Simulation timestep number
  int atomnum{};
  double box_bounds[3][2];
  double box[3][3]{};
  std::vector<Atom> atoms{};          // An array of Atom objects in a frame; initialized to be empty
  std::string coord_type{};           // "unscaled", "scaled", "unwrapped", "scaled_unwrapped"
};

struct FrameCollection {
  std::vector<Frame> frames;
  std::string coord_type;  // "unscaled", "scaled", "unwrapped", or "scaled_unwrapped"
};

// In core/trajectory.h or base.h
struct ParticleTrajectory {
  std::string type;
  std::vector<double> x, y, z;  // positions over all frames
};
