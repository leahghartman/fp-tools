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
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <unordered_map>
#include <unordered_set>

#include "core/math.h"
#include "analysis/msd.h"
#include "io/logger.h"

#include <numeric>


// Include OpenMP only if the user turns it on
#ifdef USE_OMP
#include <omp.h>
#endif

/*
* Constructor for MSDAccumulator
*
* Initializes vectors to store summed squared displacements (msd_sum) and
* the number of samples contributing to each lag (counts).
*
* Parameters:
*   - max_lag: maximum lag time (number of frames) to calculate MSD over
*/
MSDAccumulator::MSDAccumulator(int max_lag, int num_groups)
  : finalized(false) {
  msd_sum.resize(num_groups, std::vector<double>(max_lag + 1, 0.0));
  counts.resize(num_groups, std::vector<int>(max_lag + 1, 0));
}


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
void MSDAccumulator::accumulate(int group, int lag, double dr2) {
  if (group >= msd_sum.size() || lag >= msd_sum[group].size()) {
    throw std::runtime_error("MSDAccumulator::accumulate - Invalid group or lag index");
  }
  msd_sum[group][lag] += dr2;
  counts[group][lag] += 1;
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
  for (size_t g {0}; g < msd_sum.size(); ++g) {
    for (size_t i {0}; i < msd_sum[g].size(); ++i) {
      if (counts[g][i] > 0) {
        msd_sum[g][i] /= counts[g][i];
      }
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
void MSDAccumulator::write(const std::string& base_filename,
                           const std::vector<std::string>& group_labels,
                           double dt) const {

  for (size_t g = 0; g < msd_sum.size(); ++g) {
    // Generate file name: e.g., base + "_1_2_3.dat"
    std::string output_file = base_filename + "_" + group_labels[g] + ".dat";

    std::ofstream out(output_file);
    if (!out) {
      throw std::runtime_error("Failed to open file for writing MSD");
    }

    out << "# Fp-Tools MSD data file\n";
    out << "# Group: " << group_labels[g] << "\n";
    out << "# Time\t|\tMSD\n";
    out << std::scientific << std::setprecision(15);

    for (size_t lag = 0; lag < msd_sum[g].size(); ++lag) {
      double time = lag * dt;
      double avg_msd = msd_sum[g][lag];
      out << time << "\t" << avg_msd << "\n";
    }

    out.close();
  }
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
void MSDAccumulator::plot(const std::vector<std::vector<std::string>>& input_file_groups,
                          const std::vector<std::vector<std::string>>& label_groups,
                          const std::string& output_dir,
                          const std::string& format) const {
  const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";

  for (size_t i = 0; i < input_file_groups.size(); ++i) {
    std::string command = "python3 " + script_path + " msd";

    const auto& files = input_file_groups[i];
    const auto& labels = label_groups[i];

    for (size_t j = 0; j < files.size(); ++j) {
      command += " " + files[j] + " \"" + labels[j] + "\"";
    }

    command += " " + output_dir + " " + format;

    int result = std::system(command.c_str());
    if (result != 0) {
      std::cerr << "Warning: Failed to execute MSD plotting script for group " << i << ".\n";
    }
  }
}

MSDAccumulator compute_msd(const std::vector<Frame>& frames,
                           int max_lag,
                           std::vector<std::vector<std::string>> groups,
                           Logger* logger) {
  if (frames.empty() || groups.empty()) return MSDAccumulator(0);

  int max_allowed_lag = static_cast<int>(frames.size()) - 1;
  if (max_lag <= 0 || max_lag > max_allowed_lag) {
    max_lag = max_allowed_lag;
  }

  MSDAccumulator accumulator(max_lag, static_cast<int>(groups.size()));

  // Build all atom trajectories once
  std::unordered_map<int, std::vector<Atom>> atom_trajectories;
  for (const Frame& frame : frames) {
    for (const Atom& atom : frame.atoms) {
      atom_trajectories[atom.id].push_back(atom);
    }
  }

  if (groups.empty()) {
    groups.push_back({"__all__"});
  }

  std::size_t group_index = 0;
  for (const auto& group : groups) {
    if (group.size() == 1) {
      // Single-type group: optimize by skipping unordered_set
      const std::string& target_type = group[0];

      // Pre-filter trajectories of this type
      std::vector<const std::vector<Atom>*> matched_trajectories;
      matched_trajectories.reserve(atom_trajectories.size());

      for (const auto& [atom_id, traj] : atom_trajectories) {
        if (!traj.empty() && (target_type == "__all__" || traj[0].type == target_type)) {
          matched_trajectories.push_back(&traj);
        }
      }

      std::size_t done = 0;
      std::size_t total = matched_trajectories.size();

      for (int idx = 0; idx < static_cast<int>(matched_trajectories.size()); ++idx) {
        const auto& traj = *matched_trajectories[idx];
        int N = static_cast<int>(traj.size());
        for (int lag = 0; lag <= max_lag; ++lag) {
          int n_valid = N - lag;
          if (n_valid <= 0) continue;
          for (int i = 0; i < n_valid; ++i) {
            double dr2 = squared_displacement(traj[i + lag], traj[i]);
            // thread-safe accumulate here (atomic or thread-local buffers)
            accumulator.accumulate(group_index, lag, dr2);
          }
        }
        if (++done % 10 == 0 || done == total) {
          print_progress_bar(done, total, 40, "MSD (Group " + std::to_string(group_index) + ")", "particles");
        }
      }
    } else {
      std::unordered_set<std::string> type_set(group.begin(), group.end());

      // Pre-filter trajectories that match any type in the group
      std::vector<const std::vector<Atom>*> matched_trajectories;
      matched_trajectories.reserve(atom_trajectories.size());

      for (const auto& [atom_id, traj] : atom_trajectories) {
        if (!traj.empty() && type_set.count(traj[0].type)) {
          matched_trajectories.push_back(&traj);
        }
      }

      std::size_t done = 0;
      std::size_t total = matched_trajectories.size();

      for (int idx = 0; idx < static_cast<int>(matched_trajectories.size()); ++idx) {
        const auto& traj = *matched_trajectories[idx];
        int N = static_cast<int>(traj.size());
        for (int lag = 0; lag <= max_lag; ++lag) {
          int n_valid = N - lag;
          if (n_valid <= 0) continue;
          for (int i = 0; i < n_valid; ++i) {
            double dr2 = squared_displacement(traj[i + lag], traj[i]);
            accumulator.accumulate(group_index, lag, dr2);
          }
        }
        if (++done % 10 == 0 || done == total) {
          print_progress_bar(done, total, 40, "MSD (Group " + std::to_string(group_index) + ")", "particles");
        }
      }
    }
    ++group_index;
  }
  accumulator.finalize();
  return accumulator;
}
