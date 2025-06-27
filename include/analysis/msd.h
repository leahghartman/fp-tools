#pragma once
#include <map>
#include <vector>
#include <string>
#include "core/base.h"

class MSDAccumulator {
public:
  MSDAccumulator(int max_lag);   // Constructor called whenever object is instantiated

  void accumulate(int lag, double dr2);   // Adds a squared displacement value to the correct time lag
  void finalize();                        // Computes actual MSD values by dividing accumulated sums by number of contributions
  void write(const std::string& output_file) const;
  void plot(const std::string& output_file, const std::string& format) const;
    
  const std::vector<double>& get_msd() const;   // Returns a reference to the vector of MSD values
  const std::vector<int>& get_counts() const;   // Returns a reference to the number of samples (frame pairs)
  
  // Validation helpers
  static void validate_non_empty_frames(const std::vector<Frame>& frames);
  static void validate_frame_range(const std::vector<Frame>& frames, int n_start, int n_end);
  static void validate_consistent_atom_counts(const std::vector<Frame>& frames);

  // Set the path to the python plotting script (static, class-wide)
  static void set_plot_dir(const std::string& dir);

private:
  static std::string plot_dir;
  std::vector<double> msd_sum{};    // Stores the sum of squared displacements for each lag
  std::vector<int> counts{};        // Stores how many values were added for each lag (needed to compute average)
  bool finalized{false};            // Flag indicating whether finalize() has been called
};

MSDAccumulator compute_msd(const std::vector<Frame>& frames, 
                           int n_start, 
                           int n_end, 
                           int max_lag = -1);

