/**
 * -----------------------------------------------------------------------------
 * @file mfpt.cpp
 * @brief Implementation of Mean First Passage Time (MFPT) analysis tools.
 *
 * This file contains the implementation of:
 *   - FPTAccumulator methods for collecting and finalizing first-passage time 
 *     distributions.
 *   - Functions to compute FPT distributions over molecular dynamics trajectories.
 *   - Functions to write MFPT summary data and F_r(t) distributions to files.
 *   - Utilities to compute and plot the local log-derivative D(r) of MFPT data.
 *
 * Authors: Leah Hartman
 * -----------------------------------------------------------------------------
 */

#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "analysis/mfpt.h"
#include "io/logger.h"

// =============================================================================
// ===                     FPTAccumulator Class Methods                      ===
// =============================================================================

// Constructor: initialize sums and counters for first-passage times
FPTAccumulator::FPTAccumulator(int num_bins) 
  : fpt_sum(num_bins, 0.0), counts(num_bins, 0), finalized(false) {}

// Accumulate a first-passage time (FPT) into the specified lag bin
void FPTAccumulator::accumulate(int bin_index, double fpt) {
  if (bin_index >= 0 && static_cast<size_t>(bin_index) < fpt_sum.size()) {
    fpt_sum[bin_index] += fpt;
    counts[bin_index] += 1;
  }
}

// Finalize probability distribution F_r(t) by converting counts to probabilities
// using the total number of first-passage events recorded (total_counts).
void FPTAccumulator::finalize(double dt) {
  if (finalized) return; // Prevent double-finalization

  fpt_dist.resize(counts.size(), 0.0);
  total_counts = 0;
  for (int c : counts) total_counts += c; // Sum over all lag bins

  if (total_counts > 0) {
    for (size_t i{0}; i < counts.size(); ++i) {
      fpt_dist[i] = counts[i] / (total_counts * dt);  // Normalize to PDF
    }
  }
  finalized = true;
}

// Accessor methods
const std::vector<double>& FPTAccumulator::get_fpt_dist() const { return fpt_dist; }
int FPTAccumulator::get_total_counts() const { return total_counts; }

// =============================================================================
// ===                     Writing FPT and MFPT Outputs                      ===
// =============================================================================

/**
 * @brief Write FPT distributions and summarized MFPT data to files.
 *
 * For each radius bin:
 *   - Writes F_r(t) distributions to individual files (ftp_N.out).
 *   - Computes MFPT statistics t1/t0 and writes to summary file mfpt.dat.
 *
 * @param output_file Path for the MFPT summary file (mfpt.dat).
 * @param mfpt_accums Vector of accumulators containing finalized FPT data.
 * @param dr Radius bin width.
 * @param dt Time step between trajectory frames.
 */
void write_fpt(const std::string& mfpt_file,
               const std::string& fpt_dir,
               const std::vector<FPTAccumulator>& fpt_accums,
               double dr,
               double dt) {
  std::ofstream mfpt_out(mfpt_file);
  mfpt_out << "# mean first passage time T(r) vs r\n";
  mfpt_out << "# 1) r , 2) MFPT , 3) t1 , 4) t2 , 5) counts\n";

  // Loop through all radius thresholds
  for (size_t radius_idx{0}; radius_idx < fpt_accums.size(); ++radius_idx) {
    const auto& accum = fpt_accums[radius_idx];
    const auto& fpt_dist = accum.get_fpt_dist();
    int counts = accum.get_total_counts();

    double r = dr * (radius_idx + 1);  // Midpoint radius for this bin

    // Write individual F_r(t) probability distribution file for this radius
    std::string fpt_filename = fpt_dir + "/fpt_" + std::to_string(radius_idx + 1) + ".out";
    std::ofstream fpt_out(fpt_filename);
    fpt_out << "# Probability distribution F_r(t) of first passages at distance r = " << r << "\n";
    fpt_out << "# 1) t , 2) F_r(t)\n";
    fpt_out << std::scientific << std::setprecision(15);

    // Accumulate integrals for MFPT calculation
    double t0=0.0, t1=0.0, t2=0.0;
    for (size_t ibin=0; ibin<fpt_dist.size(); ++ibin) {
      double time = (ibin+1) * dt;
      double frt = fpt_dist[ibin];
      fpt_out << time << "\t" << frt << "\n";

      t0 += frt * dt;                // Zeroth moment (normalization)
      t1 += time * frt * dt;         // First moment (mean)
      t2 += time * time * frt * dt;  // Second moment
    }
    fpt_out.close();

    // Compute MFPT as <t> = t1 / t0
    double mfpt = (counts > 0 && t0 > 0) ? t1 / t0 : 0.0;

    // Write summary entry to mfpt.dat
    mfpt_out << std::scientific << std::setprecision(15)
             << r << "\t" << mfpt << "\t" << t1 << "\t" << t2 << "\t" << counts << "\n";
  }
  mfpt_out.close();
}

// =============================================================================
// ===               Compute FPT Distributions over Trajectory               ===
// =============================================================================

/**
 * @brief Accumulate FPT distributions across the trajectory frames.
 *
 * For each radius bin, accumulates first passage times for all particles,
 * starting from every possible initial frame, and records lag times at which
 * particles first cross the radius threshold.
 *
 * @param frames Trajectory frames to process.
 * @param dr Radius bin width.
 * @param r_max Maximum radius threshold.
 * @param dt Time step between frames.
 * @return Vector of MFPTAccumulator objects with accumulated FPT data.
 */
std::vector<FPTAccumulator> compute_fpt(const std::vector<Frame>& frames,
                                        double dr,
                                        double r_max,
                                        double dt) {
  size_t n_frames = frames.size();
  size_t num_atoms = frames[0].atoms.size();
  int radbins = static_cast<int>(std::ceil(r_max / dr));
  int max_lag_steps = static_cast<int>(n_frames - 1);

  // Create one accumulator per radius threshold
  std::vector<FPTAccumulator> fpt_accums;
  for (int i{0}; i < radbins; ++i) {
    fpt_accums.emplace_back(max_lag_steps + 1);
  }

  size_t total = n_frames;
  size_t done = 0;

  // Loop over every possible starting frame (time origin t0)
  for (size_t t0{0}; t0 < n_frames; ++t0) {
    const Frame& frame_start = frames[t0];

    // Record initial positions of all atoms at t0
    std::vector<double> ref_x(num_atoms), ref_y(num_atoms), ref_z(num_atoms);
    for (size_t atom_idx{0}; atom_idx < num_atoms; ++atom_idx) {
      ref_x[atom_idx] = frame_start.atoms[atom_idx].x;
      ref_y[atom_idx] = frame_start.atoms[atom_idx].y;
      ref_z[atom_idx] = frame_start.atoms[atom_idx].z;
    }

    // Track first-passage state: whether each atom has crossed each radius
    std::vector<std::vector<bool>> crossed(radbins, std::vector<bool>(num_atoms, false));

    // For each later frame after t0, check for first-passage events
    for (size_t t{t0 + 1}; t < n_frames; ++t) {
      const Frame& frame_current = frames[t];
      int lag_steps = static_cast<int>(t - t0);

      for (size_t atom_idx{0}; atom_idx < num_atoms; ++atom_idx) {
        double dx = frame_current.atoms[atom_idx].x - ref_x[atom_idx];
        double dy = frame_current.atoms[atom_idx].y - ref_y[atom_idx];
        double dz = frame_current.atoms[atom_idx].z - ref_z[atom_idx];
        double dr2 = dx * dx + dy * dy + dz * dz;

        for (int radius_idx{0}; radius_idx < radbins; ++radius_idx) {
          if (crossed[radius_idx][atom_idx]) continue;  // Already recorded

          double current_r = dr * (radius_idx + 1);
          double current_r2 = current_r * current_r;

          if (dr2 >= current_r2) {
            fpt_accums[radius_idx].accumulate(lag_steps, lag_steps);
            crossed[radius_idx][atom_idx] = true;
            break;  // Stop checking larger radii once the first is crossed
          }
        }
      }
    }
    if (++done % 10 == 0 || done == total) {
      print_progress_bar(done, total, 40, "MFPT", "origins");
    }
  }

  // Finalize the probability distribution for each radius
  for (auto& accum : fpt_accums) {
    accum.finalize(dt);
  }
  return fpt_accums;
}

// std::vector<FPTAccumulator> compute_fpt(const std::vector<Frame>& frames,
//                                         double dr,
//                                         double r_max,
//                                         double dt) {
//     if (frames.empty()) return {};
//
//     const size_t n_bins = static_cast<size_t>(r_max / dr);
//
//     // One accumulator per bin
//     std::vector<FPTAccumulator> accums;
//     accums.reserve(n_bins);
//     for (size_t i = 0; i < n_bins; ++i)
//         accums.emplace_back(n_bins);  // assumes FPTAccumulator(num_bins)
//
//     // Map atom ID to its initial position (t = 0)
//     const Frame& first_frame = frames[0];
//     std::unordered_map<std::size_t, std::tuple<double, double, double>> r0_map;
//     for (const auto& atom : first_frame.atoms)
//         r0_map[atom.id] = {atom.x, atom.y, atom.z};
//
//     // For each atom, keep track of which bins it has crossed
//     std::unordered_map<std::size_t, std::vector<bool>> crossed;
//     for (const auto& [id, _] : r0_map)
//         crossed[id] = std::vector<bool>(n_bins, false);
//
//     // Progress bar setup
//     const size_t total = (frames.size() - 1) * first_frame.atoms.size();
//     size_t processed = 0;
//
//     // Main loop over frames (excluding t=0)
//     for (size_t t = 1; t < frames.size(); ++t) {
//         double time = t * dt;
//         const auto& atoms = frames[t].atoms;
//
//         for (const auto& atom : atoms) {
//             ++processed;
//             if (processed % 10 == 0 || processed == total)
//                 print_progress_bar(processed, total, 40, "FPT accumulation");
//
//             auto r0_it = r0_map.find(atom.id);
//             if (r0_it == r0_map.end()) continue;
//
//             const auto& [x0, y0, z0] = r0_it->second;
//             double dx = atom.x - x0;
//             double dy = atom.y - y0;
//             double dz = atom.z - z0;
//             double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
//
//             size_t bin = static_cast<size_t>(dist / dr);
//             if (bin >= n_bins) continue;
//
//             // If atom hasn't yet crossed this bin
//             // if (!crossed[atom.id][bin]) {
//             //     accums[bin].accumulate(bin, time);
//             //     crossed[atom.id][bin] = true;
//             // }
//             if (dist > dr * (bin + 1)) {
//               if (!crossed[atom.id][bin]) {
//                 crossed[atom.id][bin] = true;
//                 accums[bin].accumulate(bin, time);
//               }
//             }
//         }
//     }
//
//     std::cout << std::endl;
//     return accums;
// }


// =============================================================================
// ===                         Plotting Utilities                            ===
// =============================================================================

/**
 * @brief Execute Python script to plot MFPT data.
 */
void plot_mfpt(const std::string& summary_file,
               const std::string& output_dir,
               const std::string& format) {
  namespace fs = std::filesystem;
  fs::path abs_summary_file = fs::absolute(summary_file);
  fs::path abs_output_dir = fs::absolute(output_dir);

  const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";
  std::string command = "python3 " + script_path + " mfpt " +
                        abs_summary_file.string() + " " +
                        abs_output_dir.string() + " " +
                        format;

  int result = std::system(command.c_str());
}

// =============================================================================
// ===               Log-Derivative Computation D(r) and Plot                ===
// =============================================================================

/**
 * @brief Compute local log-derivative D(r) from MFPT summary file.
 */
void compute_dofr(const std::string& mfpt_file, const std::string& output_file) {
  std::ifstream infile(mfpt_file);
  std::vector<double> r_vals, t_vals;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    double r, t;
    if (!(iss >> r >> t)) continue;
    r_vals.push_back(r);
    t_vals.push_back(t);
  }
  infile.close();

  std::ofstream outfile(output_file);
  outfile << "# Local log-derivative D(r)\n";
  outfile << "# 1) r (midpoint), 2) D(r)\n";

  for (size_t i = 1; i < r_vals.size(); ++i) {
    double r1 = r_vals[i-1], t1 = t_vals[i-1];
    double r2 = r_vals[i],   t2 = t_vals[i];

    if (r1 <= 0 || r2 <= 0 || t1 <= 0 || t2 <= 0) continue;

    double derivative = (std::log(t2) - std::log(t1)) / (std::log(r2) - std::log(r1));
    double r_mid = 0.5 * (r1 + r2);
    outfile << r_mid << " " << derivative << "\n";
  }
  outfile.close();
}

/**
 * @brief Execute Python script to plot D(r) data.
 */
void plot_dofr(const std::string& log_derivative_file,
               const std::string& output_dir,
               const std::string& format) {
  namespace fs = std::filesystem;
  fs::path abs_log_file = fs::absolute(log_derivative_file);
  fs::path abs_output_dir = fs::absolute(output_dir);

  const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";
  std::string command = "python3 " + script_path + " dlog " +
                        abs_log_file.string() + " " +
                        abs_output_dir.string() + " " +
                        format;

  int result = std::system(command.c_str());
}
