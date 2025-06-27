/*
* -----------------------------------------------------------------------------
* Mean Squared Displacement (MSD) Accumulator
* -------------------------------------------
* This file implements the MSDAccumulator class for computing Mean Squared
* Displacement (MSD) from molecular dynamics trajectories. It accumulates
* squared displacements for different time lags, averages them, writes the
* results to disk, and can call a Python script to plot the MSD curve.
*
* Authors: leahghartman
* -----------------------------------------------------------------------------
*/

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include "core/math.h"
#include "analysis/msd.h"

/*
* Constructor for MSDAccumulator
*
* Initializes vectors to store summed squared displacements (msd_sum) and
* the number of samples contributing to each lag (counts).
*
* Parameters:
*   - max_lag: maximum lag time (number of frames) to calculate MSD over
*/
MSDAccumulator::MSDAccumulator(int max_lag)
  : msd_sum(max_lag + 1, 0.0), counts(max_lag + 1, 0) {}

/*
* Accumulate a squared displacement (dr2) for a given lag index.
*
* Adds dr2 to msd_sum[lag] and increments counts[lag].
* Throws an exception if lag exceeds allocated size.
*
* Parameters:
*   - lag: index corresponding to the time lag
*   - dr2: squared displacement for this lag
*/
void MSDAccumulator::accumulate(int lag, double dr2) {
  if (lag >= msd_sum.size()) {
    throw std::runtime_error("MSDAccumulator::accumulate - Lag exceeds limit (max_lag)");
  }
  msd_sum[lag] += dr2;
  counts[lag] += 1;
}

/*
* Finalize the MSD calculation.
*
* After all accumulations, this averages the summed squared displacements
* by dividing each msd_sum[lag] by counts[lag], resulting in the mean
* squared displacement for each lag.
*
* This is called once after accumulation is complete.
*/
void MSDAccumulator::finalize() {
  if (finalized) return;
  for (size_t i{0}; i < msd_sum.size(); ++i) {
    if (counts[i] > 0) {
      msd_sum[i] /= counts[i];
    }
  }
  finalized = true;
}

/*
* Write MSD results to a text file.
*
* The output file will contain three columns per line:
*   - Lag index
*   - Mean squared displacement at that lag
*   - Number of samples (counts) contributing to that lag
*
* Parameters:
*   - output_file: path to the output file (e.g., "output/msd.dat")
*/
void MSDAccumulator::write(const std::string& output_file) const {
  std::ofstream out(output_file);
  out << "# Lag MSD Count\n";
  for (size_t i{0}; i < msd_sum.size(); ++i) {
    out << i << " " << msd_sum[i] << " " << counts[i] << "\n";
  }
  out.close();
}

// Initialize static variable to some default relative path:
std::string MSDAccumulator::plot_dir = "src/plotting/";

void MSDAccumulator::set_plot_dir(const std::string& dir) {
  plot_dir = dir;
}

/*
* Plot the MSD curve by calling an external Python script.
*
* Invokes 'plot_msd.py' with two arguments:
*   - input_file: the MSD data file (e.g., "msd.dat")
*   - format: desired output image format (e.g., "png", "pdf")
*
* Emits a warning if the Python script fails to execute.
*
* Parameters:
*   - input_file: path to the MSD data file
*   - format:     image format to generate
*/
void MSDAccumulator::plot(const std::string& input_file, const std::string& format) const {
  // Construct full script path by joining plotting_dir and script filename
  //std::string script_path = plot_dir;
  //if (!script_path.empty() && script_path.back() != '/') {
    //script_path += '/';
  //}
  const std::string script_path = std::string(SOURCE_DIR) + "/src/plotting/plot_msd.py";
  std::string command{"python3 " + script_path + " " + input_file + " " + format};
  int result{std::system(command.c_str())};

  if (result != 0) {
    std::cerr << "Warning: Failed to execute MSD plotting script.\n";
  }
}

// =========================================
// ==          Accessor Methods           ==
// =========================================

/*
* Returns the vector of mean squared displacement values after finalization.
* Declared const to prevent external modification.
*/
const std::vector<double>& MSDAccumulator::get_msd() const { return msd_sum; }

/*
* Returns the vector of counts, i.e., the number of samples contributing to
* each MSD value. Declared const to prevent external modification.
*/
const std::vector<int>& MSDAccumulator::get_counts() const { return counts; }

// =========================================
// ==         Validation Helpers          ==
// =========================================

/*
* Validates that the frame list is not empty.
*/
void MSDAccumulator::validate_non_empty_frames(const std::vector<Frame>& frames) {
  if (frames.empty()) {
    throw std::runtime_error("MSD Calculation - No frames provided.");
  }
}

/*
* Validates that n_start and n_end define a valid frame range within
* the given frames.
*/
void MSDAccumulator::validate_frame_range(const std::vector<Frame>& frames, int n_start, int n_end) {
  if (n_start < 0 || n_end > static_cast<int>(frames.size()) || n_start >= n_end - 1) {
    throw std::runtime_error("MSD Calculation - Invalid frame range");
  }
}

/*
* Validates that all frames contain the same number of atoms.
*/
void MSDAccumulator::validate_consistent_atom_counts(const std::vector<Frame>& frames) {
  int num_atoms = static_cast<int>(frames[0].atoms.size());
  for (const auto& frame : frames) {
    if (static_cast<int>(frame.atoms.size()) != num_atoms) {
      throw std::runtime_error("MSD Calculation - Inconsistent atom count across frames");
    }
  }
}

// =========================================
// ==         Main MSD Algorithm          ==
// =========================================

/*
* Compute the Mean Squared Displacement (MSD) from a series of frames.
*
* This function iterates over frame pairs separated by time lags and computes
* squared displacements for each atom, accumulating results in an MSDAccumulator.
*
* Performs basic input validation:
*   - Checks for empty frame list
*   - Ensures consistent atom counts across frames
*   - Clamps max_lag if it exceeds the available frame range
*
* Parameters:
*   - frames:   vector of simulation frames (each containing atom positions)
*   - n_start:  first frame index to include
*   - n_end:    last frame index (exclusive) to include
*   - max_lag:  maximum time lag to compute
*
* Returns:
*   - MSDAccumulator with finalized MSD values
*/
MSDAccumulator compute_msd(const std::vector<Frame>& frames, int n_start, int n_end, int max_lag) {
  MSDAccumulator::validate_non_empty_frames(frames);
  MSDAccumulator::validate_frame_range(frames, n_start, n_end);
  MSDAccumulator::validate_consistent_atom_counts(frames);
  
  int num_atoms = static_cast<int>(frames[0].atoms.size());

  // Clamp max_lag if it exceeds the allowed range
  int max_allowed_lag{n_end - n_start - 1};
  if (max_lag <= 0 || max_lag > max_allowed_lag) {
    max_lag = max_allowed_lag;
  }
  
  MSDAccumulator accumulator{max_lag};

  // Loop over lag times from 1 up to max_lag
  for (int lag{1}; lag <= max_lag; ++lag) {
    double sum_dr2{0.0};
    double count{0};

    // Iterate over all valid frame pairs separated by 'lag'
    for (int t{n_start}; t + lag < n_end; ++t) {
      const Frame& f0{frames[t]};
      const Frame& f1{frames[t + lag]};

      // Compute squared displacement for each atom
      for (int i{0}; i < num_atoms; ++i) {
        double dr2 = squared_displacement(f1.atoms[i], f0.atoms[i]);
        sum_dr2 += dr2;
        count++;
      }
    }
    if (count > 0) {
      accumulator.accumulate(lag, sum_dr2);
    }
  }
  accumulator.finalize();
  return accumulator;
}
