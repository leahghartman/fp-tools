/**
 * -----------------------------------------------------------------------------
 * @file corr.h
 * @brief Correlation function analysis using FFTW for velocity, force, etc.
 *
 * This header defines the CorrAccumulator class and helper functions for
 * computing time correlation functions such as VACF or custom cross-correlations
 * between time series extracted from simulation outputs. FFT-based convolution 
 * is used for performance.
 *
 * Supported features:
 *  - VACF (Velocity Autocorrelation Function)
 *  - Custom two-signal correlations (e.g., density-potential)
 *  - Mean subtraction and normalization
 *  - Zero-padding options (next_pow2, double, none)
 *
 * @author Leah Hartman (leahghartman) (LANL, University of Michigan)
 * -----------------------------------------------------------------------------
 */

#pragma once

#include <map>
#include <string>
#include <vector>
#include <fftw3.h>

#include "io/config_par.h"
#include "core/base.h"

/**
 * @class CorrAccumulator
 * @brief Computes built-in and custom correlation functions via FFTW.
 */
class CorrAccumulator {
public:
  // Constructor
  CorrAccumulator();
  ~CorrAccumulator();

  void add_builtin(const std::string& name, const std::string& input_dir);
  void add_custom(const std::string& name,
                  const std::vector<std::string>& files,
                  const std::vector<int>& columns,
                  const std::string& input_dir);

  void accumulate(const std::string& zero_pad = "next_pow2",
                  bool subtract_mean = true,
                  bool normalize = true);

  // Optional: write g(t) data to files
  void write(const std::string& output_dir) const;

  // Optional: plot using Python
  void plot(const std::string& output_dir, const std::string& format) const;

  std::unordered_map<std::string, std::vector<double>> results;

  void store_result(const std::string& name, const std::vector<double>& data) {
    results[name] = data;
  }
  void merge(const CorrAccumulator& other);

private:
  // Correlation computation
  std::vector<double> correlate(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::string& pad_method,
                                bool subtract_mean,
                                bool normalize);

  // FFT utilities
  size_t get_padded_size(size_t n, const std::string& method) const;

  // I/O utilities
  std::vector<double> read_column(const std::string& filename, int column_index);
  void save_result(const std::string& name, const std::vector<double>& g_t);

  // Results
  std::vector<std::string> result_names;
  std::vector<std::vector<double>> g_t_list;

  // FFT workspace (reused)
  fftw_complex *fft_a = nullptr;
  fftw_complex *fft_b = nullptr;
  fftw_complex *ifft_result = nullptr;
  fftw_plan plan_a;
  fftw_plan plan_b;
  fftw_plan plan_ifft;

  // FFTW buffer size
  size_t fft_size = 0;

  // Cleanup
  void cleanup_fftw();
};

CorrAccumulator compute_correlations(const CorrConfig& cfg,
                                     const std::vector<Frame>& frames);

// void compute_correlations(const std::string& input_dir,
//                const std::vector<std::string>& builtins,
//                const std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<int>>>& custom,
//                const std::string& zero_pad = "next_pow2",
//                bool subtract_mean = true,
//                bool normalize = true,
//                const std::string& plot_format = "png");
//
