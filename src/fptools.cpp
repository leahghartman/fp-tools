/**
 * -----------------------------------------------------------------------------
 * @file fptools.cpp
 * @brief Main driver for FP-Tools, a molecular dynamics analysis suite.
 *
 * Responsibilities:
 *   - Parse input TOML configuration
 *   - Detect input trajectory format and instantiate parser
 *   - Load frames from supported formats (LAMMPS, XYZ, VASP)
 *   - Perform configured analyses: MSD, RDF, MFPT, etc.
 *   - Log runtime metadata, generate output files, and optionally plot results
 *
 * @author Leah Hartman (leahghartman) (LANL, University of Michigan)
 * -----------------------------------------------------------------------------
 */

#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <unistd.h>
#include <iostream>
#include <filesystem>
#include <string_view>
#include <chrono>

// Core data structures and configuration parsing
#include "core/traj_prov.h"
#include "io/config_par.h"

// Input trajectory parsers
#include "io/lamptraj_par.h"
//#include "io/vasptraj_par.h"
#include "io/xyztraj_par.h"

// Analysis modules
#include "analysis/msd.h"
#include "analysis/rdf.h"
#include "analysis/mfpt.h"
#include "analysis/corr.h"

// Configuration parser and serializer for C++
#include "toml.hpp"

// Log files
#include "io/logger.h"

namespace fs = std::filesystem;
char hostname[1024];

std::string get_timestamp_string() {
  std::time_t t = std::time(nullptr);
  char buf[100];
  std::strftime(buf, sizeof(buf), "%F %T", std::localtime(&t));
  return std::string(buf);
}

std::string infer_format(const std::string& filename) {
  if (filename.find("XDATCAR") != std::string::npos) return "vasp";
  if (filename.ends_with(".xyz")) return "xyz";
  if (filename.ends_with(".lammpstrj") || filename.ends_with(".dump")) return "lammps";
  throw std::runtime_error(
    "Unrecognized file format: " + filename +
    "\nSupported formats are: XDATCAR (vasp), .xyz, .lammpstrj, .dump"
  );
}

int main(int argc, char** argv) {
  // Check for required configuration call format
  if (argc < 2) {
    std::cerr << "Usage: fptools <config.yaml>\n";
    return 1;
  }

  double full_elapsed = 0.0;
  auto full_start = std::chrono::steady_clock::now();

  std::string config_file{argv[1]};
  Config config = parse_config(config_file);

  std::string output_dir = config.output.path;
  std::string msd_dir = output_dir + "/msd";
  std::string rdf_dir = output_dir + "/rdf";
  std::string mfpt_dir = output_dir + "/mfpt";
  std::string corr_dir = output_dir + "/corr";

  int status_interval = config.output.status_interval;

  // Create output directories silently, error if fail
  if (!fs::exists(output_dir) && !fs::create_directories(output_dir)) {
    std::cerr << "Error: Could not create output directory: " << output_dir << "\n";
    return 1;
  }

  // Create a log file
  Logger logger(output_dir + "/run.log");
  gethostname(hostname, 1024);

  logger.banner("0.1.0");

  logger.sys_info(hostname, std::filesystem::current_path(), get_timestamp_string());
  logger.newline();

  // Parsing config
  logger.arrow("Parsing configuration file: " + config_file);

  // Infer trajectory format
  std::string format = infer_format(config.input.file);
  logger.arrow("Detected input trajectory format: " + format);
  logger.newline();

  // Instantiate trajectory parser
  logger.section("Trajectory Information");
  logger.newline();
  logger.arrow("Opening and parsing trajectory file: " + config.input.file);

  // ---------------------------------------------------
  // ---  Instantiate appropriate trajectory parser  ---
  // ---------------------------------------------------
  std::unique_ptr<TrajectoryProvider> traj{};

  if (format == "lammps") {
    traj = std::make_unique<LAMMPSTrajectory>(config.input.file);
  } else if (format == "xyz") {
    traj = std::make_unique<XYZTrajectory>(config.input.file);
  } else if (format == "vasp") {
    return 1;//traj = std::make_unique<VASPTrajectory>(config.input.file);
  }
  
  // -----------------------------------------------------
  // ---  Load frames according to config information  ---
  // -----------------------------------------------------
  std::vector<Frame> frames{};
  int frame_count{0};
  Frame frame{};
  std::size_t loaded_frames{0};

  int frame_span = (config.input.end_frame >= 0)
                   ? config.input.end_frame - config.input.start_frame + 1
                   : std::numeric_limits<int>::max();  // fallback if no end frame

  int expected_frames = (frame_span == std::numeric_limits<int>::max())
                        ? -1 : (frame_span + config.input.frame_interval - 1) / config.input.frame_interval;

  while (traj->next_frame(frame)) {
    if (config.input.end_frame >= 0 && frame_count > config.input.end_frame)
      break;

    if (frame_count >= config.input.start_frame &&
        ((frame_count - config.input.start_frame) % config.input.frame_interval == 0)) {
      frames.push_back(frame);
      ++loaded_frames;

      if (expected_frames > 0 &&
        (loaded_frames % 10 == 0 || loaded_frames == expected_frames)) {
        print_progress_bar(loaded_frames, expected_frames, 40, "Loading frames", "frames");
      }
    }
    ++frame_count;
  }

  double box_volume = 1.0;
  if (!frames.empty()) {
    const auto& bounds = frames[0].box_bounds;
    double lx = bounds[0][1] - bounds[0][0];
    double ly = bounds[1][1] - bounds[1][0];
    double lz = bounds[2][1] - bounds[2][0];
    box_volume = lx * ly * lz;
  }

  
  //logger.newline();
  logger.arrow("Trajectory parsing completed successfully!");
  logger.newline();
  logger.traj_info(config.input, frames);

  // ----------------------------------------
  // ---  Prepare for output and logging  ---
  // ----------------------------------------
  int verbosity{config.output.verbosity};

  // ------------------------------------
  // ---  Perform requested analyses  ---
  // ------------------------------------
  
  bool ran_msd = false;
  double msd_dt = 0.0;
  int msd_max_lag = 0;
  int msd_num_groups = 0;
  int msd_num_frames = static_cast<int>(frames.size());
  double msd_elapsed = 0.0;
  std::vector<std::vector<std::string>> msd_groups;

  if (config.analysis.msd.enabled) {
    ran_msd = true;
    msd_dt = config.properties.dt;
    msd_max_lag = config.analysis.msd.max_lag;
    msd_groups = config.analysis.msd.groups;
    msd_num_groups = static_cast<int>(msd_groups.size());
    double dt               = config.properties.dt;
    int max_lag             = config.analysis.msd.max_lag;
    std::string plot_format = config.analysis.msd.plot_format;
    const auto& groups      = config.analysis.msd.groups;
    const auto& curve_sets  = config.analysis.msd.curves;

    fs::create_directories(msd_dir);

    logger.msd_info(dt, max_lag, groups);
    auto msd_start = std::chrono::steady_clock::now();

    MSDAccumulator msd_accum(max_lag);
    msd_accum = compute_msd(frames, max_lag, groups, &logger);

    // Build group labels like "1_2_3"
    std::vector<std::string> group_labels;
    for (const auto& group : groups) {
      std::ostringstream oss;
      for (size_t i = 0; i < group.size(); ++i) {
        if (i > 0) oss << "_";
        oss << group[i];
      }
      group_labels.push_back(oss.str());
    }

    // Write .dat files with label-based names
    logger.arrow("Writing MSD data to " + msd_dir + "...");
    msd_accum.write(msd_dir + "/msd", group_labels, config.properties.dt);

    std::vector<std::vector<std::string>> input_file_groups;
    std::vector<std::vector<std::string>> label_groups;

    for (const auto& curve_group : curve_sets) {
      std::vector<std::string> file_group;
      std::vector<std::string> label_group;

      for (const auto& curve : curve_group) {
        std::string types_label;
        for (size_t i = 0; i < curve.types.size(); ++i) {
          if (i > 0) types_label += "_";
          types_label += curve.types[i];
        }
        std::string filename = msd_dir + "/msd_" + types_label + ".dat";
        file_group.push_back(filename);
        label_group.push_back(curve.label);
      }
      input_file_groups.push_back(file_group);
      label_groups.push_back(label_group);
    }

    // Plot results
    if (!plot_format.empty()) {
      logger.arrow("Generating MSD plot(s) (" + plot_format + ")...");
      msd_accum.plot(input_file_groups, label_groups, msd_dir, plot_format);
      logger.arrow("Plot(s) successfully saved to " + msd_dir + "/!");
    }
    auto msd_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = msd_end - msd_start;
    msd_elapsed = elapsed_seconds.count();
    logger.newline();
  }


  // Radial Density Function (RDF) calculation
  std::vector<std::vector<std::string>> rdf_pairs = config.analysis.rdf.pairs;
  int rdf_num_frames = static_cast<int>(frames.size());
  double rdf_dt = config.properties.dt;
  double rdf_elapsed = 0.0;
  double rdf_dr = config.analysis.rdf.dr;
  int rdf_num_bins = config.analysis.rdf.num_bins;
  double rdf_r_max = config.analysis.rdf.r_max;
  bool ran_rdf = false;

  if (config.analysis.rdf.enabled) {
    ran_rdf = true;
    auto rdf_start = std::chrono::steady_clock::now(); // ⏱️ Start timing

    int num_bins                                {config.analysis.rdf.num_bins};
    double dr                                   {config.analysis.rdf.dr};
    double r_max                                {config.analysis.rdf.r_max};
    std::vector<std::vector<std::string>> pairs {config.analysis.rdf.pairs};
    std::string plot_format                     {config.analysis.rdf.plot_format};

    fs::create_directories(rdf_dir);

    logger.rdf_info(dr, num_bins, r_max, pairs);

    // --- Ensure atom types exist, or default to "1"
    bool has_missing_types = frames[0].atoms.empty() || frames[0].atoms[0].type.empty();
    if (has_missing_types) {
      for (auto& frame : const_cast<std::vector<Frame>&>(frames)) {
        for (auto& atom : frame.atoms) {
          atom.type = "1";
        }
      }
    }

    // --- Instantiate and compute RDF
    RDFAccumulator rdf_accum = (num_bins > 0) ? RDFAccumulator(static_cast<size_t>(num_bins), r_max)
                                            : RDFAccumulator(static_cast<double>(dr), r_max);

    if (num_bins > 0) {
      rdf_accum = compute_rdf(frames, static_cast<std::size_t>(num_bins), r_max, box_volume, pairs, &logger);
    } else if (dr > 0.0) {
      rdf_accum = compute_rdf(frames, dr, r_max, box_volume, pairs, &logger);
    } else {
      throw std::runtime_error("Neither num_bins nor dr is properly specified in config.");
    }

    logger.arrow("Writing RDF data to " + rdf_dir + "...");
    rdf_accum.write(rdf_dir);

    // --- Prepare label_groups and input_file_groups
    std::vector<std::vector<std::string>> label_groups;
    std::vector<std::vector<std::string>> input_file_groups;

    for (const auto& pair : pairs) {
      std::ostringstream oss;
      for (size_t i = 0; i < pair.size(); ++i) {
        if (i > 0) oss << "_";
        oss << pair[i];
      }
      std::string label = oss.str();
      std::string filename = rdf_dir + "/rdf-" + label + ".dat";

      label_groups.push_back({label});
      input_file_groups.push_back({filename});
    }

    // --- Plot RDF if enabled
    if (!plot_format.empty()) {
      logger.arrow("Generating RDF plot(s) (" + plot_format + ")...");
      rdf_accum.plot(input_file_groups, label_groups, rdf_dir, plot_format);
      logger.arrow("Plot(s) successfully saved to " + rdf_dir + "/!");
    }

    auto rdf_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = rdf_end - rdf_start;
    rdf_elapsed = elapsed_seconds.count();
    logger.newline();
  }
  
  bool ran_mfpt = false;
  double mfpt_elapsed = 0.0;
  int mfpt_num_radii = 0;
  int mfpt_num_frames = static_cast<int>(frames.size());
  double mfpt_dt = config.analysis.mfpt.dt;
  double mfpt_dr = config.analysis.mfpt.dr;
  double mfpt_r_max = config.analysis.mfpt.r_max;

  if (config.analysis.mfpt.enabled) {
    ran_mfpt = true;

    // Extract parameters
    double dr               = config.analysis.mfpt.dr;
    double r_max            = config.analysis.mfpt.r_max;
    double dt               = config.analysis.mfpt.dt;
    std::string plot_format = config.analysis.mfpt.plot_format;

    // Output paths
    std::string fpt_dir     = mfpt_dir + "/fpt";
    std::string mfpt_file   = mfpt_dir + "/mfpt.dat";
    std::string dr_file     = output_dir + "/D_r.dat";

    fs::create_directories(fpt_dir);

    // Info + timer
    logger.mfpt_info(dr, r_max, dt);
    auto mfpt_start = std::chrono::steady_clock::now();

    // Compute and write MFPT data
    std::vector<FPTAccumulator> fpt_accums = compute_fpt(frames, dr, r_max, dt);
    logger.arrow("Writing MFPT data to " + mfpt_dir + "...");
    write_fpt(mfpt_file, fpt_dir, fpt_accums, dr, dt, box_volume, frames[0].atomnum);

    // Compute D(r)
    logger.arrow("Computing D(r) from MFPT...");
    compute_dofr(mfpt_file, dr_file);

    logger.arrow("Plotting D(r)...");
    plot_dofr(dr_file, mfpt_dir, "svg");

    if (!plot_format.empty()) {
      logger.arrow("Generating MFPT plot(s) (" + plot_format + ")...");
      plot_mfpt(mfpt_file, mfpt_dir, plot_format);
      logger.arrow("Plot(s) successfully saved to " + mfpt_dir + "/!");
    }

    auto mfpt_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = mfpt_end - mfpt_start;
    double mfpt_elapsed = elapsed_seconds.count();
    mfpt_num_radii = static_cast<int>(mfpt_r_max / mfpt_dr);

    logger.newline();
  }

  if (config.analysis.corr.enabled) {
    bool normalize {config.analysis.corr.normalize};
    bool subtract_mean {config.analysis.corr.subtract_mean};
    std::string zero_pad {config.analysis.corr.zero_pad};
    std::string plot_format {config.analysis.corr.plot_format};
    std::vector<std::string> builtins {config.analysis.corr.builtins};
    std::vector<CustomCorr> custom {config.analysis.corr.custom};

    std::vector<std::string> custom_names;
    for (const auto& c : custom) {
      custom_names.push_back(c.name);
    }

    fs::create_directories(corr_dir);
    std::string corr_file{corr_dir + "/corr.dat"};

    logger.corr_info(builtins, custom_names, zero_pad, normalize, subtract_mean);

    CorrAccumulator corr_accum = compute_correlations(config.analysis.corr, frames);

    logger.arrow("Writing correlation data to " + corr_file);
    corr_accum.write(corr_dir);

    // Optionally plot correlation results using configured plot format
    if (!plot_format.empty()) {
      logger.arrow("Generating correlation plot (" + plot_format + ")...");
      //corr_accum.plot(corr_file, corr_dir, plot_format);
    logger.arrow("Plot successfully saved to " + corr_dir + "/!");
    }
    logger.newline();
  }
  
  auto full_end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = full_end - full_start;

  logger.section("Summary");
  logger.newline();
  if (ran_msd) {
    logger.msd_summ(msd_groups, msd_num_frames, msd_dt, msd_elapsed);
  }
  if (ran_rdf) {
    logger.rdf_summ(rdf_pairs,
                rdf_num_frames,
                rdf_dt,
                rdf_elapsed,
                rdf_dr,
                rdf_num_bins,
                rdf_r_max);
  }
  if (ran_mfpt) {
    logger.mfpt_summ(mfpt_num_radii, mfpt_num_frames, mfpt_dr, mfpt_r_max, mfpt_dt, mfpt_elapsed);
  }


  std::ostringstream oss;
  oss << "[ TOTAL ] Total runtime: " << std::fixed << std::setprecision(2)
      << elapsed_seconds.count() << " seconds";
  logger.note(oss.str());
  logger.newline();
  return 0;
}

