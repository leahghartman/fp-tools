
#include "analysis/corr.h"
#include "core/base.h"
#include "io/config_par.h"
#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

CorrAccumulator::CorrAccumulator() = default;
CorrAccumulator::~CorrAccumulator() = default;

// Helper: subtract mean in-place
void subtract_mean(std::vector<double> &data) {
  double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
  for (auto &val : data)
    val -= mean;
}

// Helper: zero pad to nearest power of 2
size_t next_pow2(size_t n) {
  size_t pow2 = 1;
  while (pow2 < n)
    pow2 <<= 1;
  return pow2;
}

void CorrAccumulator::merge(const CorrAccumulator &other) {
  for (const auto &[key, val] : other.results) {
    results[key] = val;
  }
}

// Helper: read a specific column from file (1-based indexing)
std::vector<double> read_column(const std::string &path, size_t col_index) {
  std::ifstream infile(path);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file: " + path);
  }

  std::vector<double> column_data;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    double val;
    size_t current_col = 1;
    while (iss >> val) {
      if (current_col == col_index) {
        column_data.push_back(val);
        break;
      }
      ++current_col;
    }
  }
  return column_data;
}

// Core FFT-based correlation
std::vector<double> correlate(const std::vector<double> &a,
                              const std::vector<double> &b,
                              const CorrConfig &cfg) {
  size_t N = a.size();
  size_t M = b.size();

  size_t pad_size = cfg.zero_pad == "next_pow2" ? next_pow2(2 * N) : 2 * N;
  std::vector<double> A(pad_size, 0.0);
  std::vector<double> B(pad_size, 0.0);
  std::copy(a.begin(), a.end(), A.begin());
  std::copy(b.begin(), b.end(), B.begin());

  fftw_complex *fft_a =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (pad_size / 2 + 1));
  fftw_complex *fft_b =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (pad_size / 2 + 1));
  fftw_complex *fft_out =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (pad_size / 2 + 1));

  fftw_plan plan_a =
      fftw_plan_dft_r2c_1d(pad_size, A.data(), fft_a, FFTW_ESTIMATE);
  fftw_plan plan_b =
      fftw_plan_dft_r2c_1d(pad_size, B.data(), fft_b, FFTW_ESTIMATE);
  fftw_plan plan_back = fftw_plan_dft_c2r_1d(pad_size, fft_out, A.data(),
                                             FFTW_ESTIMATE); // reuse A

  fftw_execute(plan_a);
  fftw_execute(plan_b);

  for (size_t i = 0; i < pad_size / 2 + 1; ++i) {
    double re = fft_a[i][0] * fft_b[i][0] + fft_a[i][1] * fft_b[i][1];
    double im = fft_a[i][1] * fft_b[i][0] - fft_a[i][0] * fft_b[i][1];
    fft_out[i][0] = re;
    fft_out[i][1] = im;
  }

  fftw_execute(plan_back);

  std::vector<double> result(N);
  for (size_t i = 0; i < N; ++i) {
    double val = A[i] / pad_size;
    result[i] = cfg.normalize ? val / result[0] : val;
  }

  fftw_destroy_plan(plan_a);
  fftw_destroy_plan(plan_b);
  fftw_destroy_plan(plan_back);
  fftw_free(fft_a);
  fftw_free(fft_b);
  fftw_free(fft_out);

  return result;
}

CorrAccumulator compute_builtin(const CorrConfig &cfg,
                                const std::vector<Frame> &frames) {
  CorrAccumulator acc;
  size_t num_atoms = frames[0].atoms.size();
  size_t num_frames = frames.size();

  if (std::find(cfg.builtins.begin(), cfg.builtins.end(), "vacf") != cfg.builtins.end()) {
    std::vector<double> vacf(num_frames, 0.0);

    // Temporary vector to accumulate results from each atom
    std::vector<double> vacf_accum(num_frames, 0.0);

    for (size_t i = 0; i < num_atoms; ++i) {
      // Extract velocity components for atom i
      std::vector<double> vx(num_frames), vy(num_frames), vz(num_frames);
      for (size_t t = 0; t < num_frames; ++t) {
        vx[t] = frames[t].atoms[i].vx;
        vy[t] = frames[t].atoms[i].vy;
        vz[t] = frames[t].atoms[i].vz;
      }

      // Compute autocorrelations per component
      std::vector<double> c_x = correlate(vx, vx, cfg);
      std::vector<double> c_y = correlate(vy, vy, cfg);
      std::vector<double> c_z = correlate(vz, vz, cfg);

      // Sum components to get total VACF for this atom
      for (size_t lag = 0; lag < num_frames; ++lag) {
        vacf_accum[lag] += c_x[lag] + c_y[lag] + c_z[lag];
      }
    }

    // Average over all atoms
    for (size_t lag = 0; lag < num_frames; ++lag) {
      vacf[lag] = vacf_accum[lag] / static_cast<double>(num_atoms);
    }

    acc.store_result("vacf", vacf);
  }
  return acc;
}

CorrAccumulator compute_correlations(const CorrConfig &cfg,
                                     const std::vector<Frame> &frames) {
  CorrAccumulator acc;

  // --- Built-in functions (vacf, vcf) ---
  if (!cfg.builtins.empty()) {
    CorrAccumulator builtin_acc = compute_builtin(cfg, frames);
    acc.merge(builtin_acc); // You'll implement merge()
  }

  // --- Custom file/column pair correlations ---
  for (const auto &pair : cfg.custom) {
    std::vector<double> a, b;

    if (pair.files.size() == 1) {
      // Autocorrelation: use same file and column for both
      a = read_column(pair.files[0], pair.columns[0]);
      b = a;
    } else {
      a = read_column(pair.files[0], pair.columns[0]);
      b = read_column(pair.files[1], pair.columns[1]);
    }

    if (cfg.subtract_mean) {
      subtract_mean(a);
      // Only subtract mean from b if it's not a reference to a
      if (&a != &b)
        subtract_mean(b);
    }
    auto result = correlate(a, b, cfg);
    acc.store_result(pair.name, result);
  }
  return acc;
}

void CorrAccumulator::write(const std::string &output_dir) const {
  namespace fs = std::filesystem;
  fs::create_directories(output_dir);

  for (const auto &[name, data] : results) {
    std::string outpath = output_dir + "/" + name + "corr.dat";
    std::ofstream out(outpath);
    for (size_t i = 0; i < data.size(); ++i) {
      out << i << " " << data[i] << "\n";
    }
  }
}

void CorrAccumulator::plot(const std::string &output_dir,
                           const std::string &format) const {
  namespace fs = std::filesystem;
  fs::path abs_output_dir = fs::absolute(output_dir);
  const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";

  std::string command = "python3 " + script_path + " corr " +
                        abs_output_dir.string() + " " + format;

  int result = std::system(command.c_str());
  if (result != 0) {
    std::cerr << "Warning: Failed to run correlation plotting script.\n";
  }
}
