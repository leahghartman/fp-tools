/**
 * -----------------------------------------------------------------------------
 * @file rdf.h
 * @brief Radial Distribution Function (RDF) Accumulator and Analysis
 *
 * This header defines the RDFAccumulator class and helper functions for 
 * computing and analyzing radial distribution functions g(r) from molecular 
 * dynamics trajectories. It supports histogram accumulation for arbitrary 
 * type pairs, normalization, output, and optional plotting via an external 
 * Python script.
 *
 * Provided functions include:
 *  - compute_rdf: process a trajectory and return RDFAccumulator object
 *  - write: write g(r) results to file
 *  - plot: call Python script to generate RDF plots
 *
 * @author Leah Hartman (leahghartman) (LANL, University of Michigan)
 * -----------------------------------------------------------------------------
 */

#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include "io/logger.h"
#include "core/base.h"

/**
 * @class RDFAccumulator
 * @brief Collects and normalizes pair-separation histograms into g(r).
 */
class RDFAccumulator {
public:
  // Constructors
  RDFAccumulator(double dr, double r_max);
  RDFAccumulator(size_t num_bins, double r_max);

  // Main accumulation function; user provides pairs specifying type combinations
  void accumulate(const std::vector<Frame>& frames,
                  const std::vector<std::vector<std::string>>& pairs, Logger* logger);

  // Normalize accumulated histograms to produce g(r)
  void normalize(double box_volume, size_t num_frames,
                 const std::unordered_map<std::string, size_t>& atoms_per_type_counts);

  // Write RDF data files for each pair
  void write(const std::string& output_dir /* std::vector<std::vector<std::string>> pairs */) const;

  // Optional: plot RDF curves (calls external python script)
  void plot(const std::vector<std::vector<std::string>>& input_file_groups,
                          const std::vector<std::vector<std::string>>& label_groups,
                          const std::string& output_dir,
                          const std::string& format) const;

  void set_include_intramolecular_pairs(bool include) {
    include_intramolecular_pairs_ = include;
  }

  // Accessors for r bins and g(r) curves
  const std::vector<double>& get_r() const { return r_vals; }
  const std::vector<std::vector<double>>& get_gr_per_pair() const { return g_r; }

private:
  // Private data members
  double dr{};
  size_t num_bins{};
  double r_max{};
  double box_volume{};  // set during normalization or externally if needed

  std::vector<double> r_vals;

  // One histogram and g_r vector per pair
  std::vector<std::vector<double>> histogram;
  std::vector<std::vector<double>> g_r;

  bool include_intramolecular_pairs_ = false;

  int nanoflann_leaf_max_size = 10; // Default, can be tuned


  // The pairs of types for which RDF is calculated
  std::vector<std::vector<std::string>> pairs;

  // Private helper functions
  void clear_histograms();
  void increment_histogram(size_t pair_index, size_t bin, double weight);
};

// Factory functions to compute RDF by specifying dr or num_bins
RDFAccumulator compute_rdf(const std::vector<Frame>& frames,
                           double dr,
                           double r_max,
                           double box_volume, 
                           std::vector<std::vector<std::string>> pairs,
                           Logger* logger);

RDFAccumulator compute_rdf(const std::vector<Frame>& frames,
                           size_t num_bins,
                           double r_max,
                           double box_volume, 
                           std::vector<std::vector<std::string>> pairs,
                           Logger* logger);
