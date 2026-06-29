// ==============================================================================
// demo_trex_gvs_06_mc_sim_ar1_blocks.cpp
// ==============================================================================
/**
 * @file demo_trex_gvs_06_mc_sim_ar1_blocks.cpp
 *
 * @brief MC simulation: T-Rex+GVS on AR(1)-within-block DGP.
 *
 * @details
 * This benchmark uses four grouped blocks of sizes {20, 50, 80, 65} with
 * AR(1) covariance structure within each block. Blocks 0, 1, and 2 are active
 * with beta = 3, while block 3 is an inactive AR(1) trap block.
 *
 * Three parts:
 * 1. SNR sweep at fixed rho = 0.85
 * 2. rho sweep at fixed SNR = 2.0
 * 3. 2-D SNR × rho sweep
 *
 * DGP:
 * - make_ar1_blocks_dgp()
 * - within-block covariance Sigma(i,j) = rho^|i-j|
 * - block sizes = {20, 50, 80, 65}
 * - active blocks = {0, 1, 2}
 * - inactive trap block = {3}
 * - support size s = 150
 *
 * Shared selector settings:
 * - K = 20
 * - tFDR = 0.1
 * - corr_max = 0.98
 * - single-linkage HAC
 * - lambda2 = CV_1SE_CCD
 *
 * Methods compared:
 * - EN (TENET)
 * - EN (TENET_AUG)
 * - IEN
 *
 * MC:
 * - 200 trials per grid point
 */
// ==============================================================================

#include "trex_gvs_dgps.hpp"
#include "trex_gvs_simulation_utils.hpp"

// ==============================================================================

using namespace gvs_demo;

// =============================================================================
// Preset and method bundle
// =============================================================================

struct AR1BlocksPreset {
  std::string scenario_tag;
  std::string file_stem_prefix;

  double fixed_rho_for_snr = 0.85;
  double fixed_snr_for_rho = 2.0;

  std::vector<double> snr_grid_1d;
  std::vector<double> rho_grid_1d;
  std::vector<double> snr_grid_2d;
  std::vector<double> rho_grid_2d;
};

struct AR1BlocksMethodSet {
  TRexControlParameter trex_ctrl;
  TRexGVSControlParameter en;
  TRexGVSControlParameter en_aug;
  TRexGVSControlParameter ien;
};

// =============================================================================
// Preset builder
// =============================================================================

static AR1BlocksPreset make_ar1_blocks_preset() {
  AR1BlocksPreset p;
  p.scenario_tag = "AR1-Blocks";
  p.file_stem_prefix = "gvs_ar1_blocks";

  p.fixed_rho_for_snr = 0.85;
  p.fixed_snr_for_rho = 2.0;

  p.snr_grid_1d = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_1d = {0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                   0.90, 0.95, 0.99};
  p.snr_grid_2d = {0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_2d = {0.30, 0.50, 0.70, 0.85, 0.90, 0.99};

  return p;
}

// =============================================================================
// Shared method controls
// =============================================================================

static AR1BlocksMethodSet make_method_set(const GVSSimConfig& cfg) {
  AR1BlocksMethodSet ms;
  ms.trex_ctrl = make_gvs_trex_control(cfg.K);

  ms.en.gvs_type = GVSType::EN;
  ms.en.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
  ms.en.corr_max = cfg.corr_max;
  ms.en.hc_linkage = hac::LinkageMethod::Single;
  ms.en.en_solver = ENSolverType::TENET;

  ms.en_aug = ms.en;
  ms.en_aug.en_solver = ENSolverType::TENET_AUG;

  ms.ien.gvs_type = GVSType::IEN;
  ms.ien.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
  ms.ien.corr_max = cfg.corr_max;
  ms.ien.hc_linkage = hac::LinkageMethod::Single;

  return ms;
}

// =============================================================================
// Part 1 — SNR sweep
// =============================================================================

/**
 * @brief Run an MC SNR sweep on the AR(1)-blocks DGP at fixed rho = 0.85.
 */
static void run_ar1_part1_snr_sweep(const GVSSimConfig& cfg,
                                    const AR1BlocksPreset& preset,
                                    const AR1BlocksMethodSet& ms) {
  const double fixed_rho = preset.fixed_rho_for_snr;

  GVSDGPFactory snr_fn = [&cfg, fixed_rho](double snr, unsigned seed) {
    return make_ar1_blocks_dgp(cfg.n, cfg.p, snr, fixed_rho, seed);
  };

  cdiag::print_section_header(
      "Part 1: SNR Sweep | " + preset.scenario_tag +
      "\nn=" + std::to_string(cfg.n) + " p=" + std::to_string(cfg.p) +
      " s=150 MC=" + std::to_string(cfg.num_MC) +
      " tFDR=" + std::to_string(cfg.tFDR).substr(0, 3) +
      " rho=0.85 (fixed)");

  auto en = run_gvs_snr_sweep(
      snr_fn, preset.snr_grid_1d, cfg, ms.en, ms.trex_ctrl, "EN");
  auto en_aug = run_gvs_snr_sweep(
      snr_fn, preset.snr_grid_1d, cfg, ms.en_aug, ms.trex_ctrl, "EN+AUG");
  auto ien = run_gvs_snr_sweep(
      snr_fn, preset.snr_grid_1d, cfg, ms.ien, ms.trex_ctrl, "IEN");

  print_mc_snr_table(
      preset.scenario_tag,
      preset.snr_grid_1d,
      en,
      en_aug,
      ien,
      cfg,
      preset.file_stem_prefix + "_snr");
}

// =============================================================================
// Part 2 — rho sweep
// =============================================================================

/**
 * @brief Run an MC rho sweep on the AR(1)-blocks DGP at fixed SNR.
 */
static void run_ar1_part2_rho_sweep(const GVSSimConfig& cfg,
                                    const AR1BlocksPreset& preset,
                                    const AR1BlocksMethodSet& ms) {
  GVSRhoDGPFactory rho_fn = [&cfg](double rho, unsigned seed) {
    return make_ar1_blocks_dgp(cfg.n, cfg.p, cfg.snr, rho, seed);
  };

  cdiag::print_section_header(
      "Part 2: rho Sweep | " + preset.scenario_tag +
      " (SNR=" + std::to_string(cfg.snr).substr(0, 3) + ")");

  auto en = run_gvs_rho_sweep(
      rho_fn, preset.rho_grid_1d, cfg, ms.en, ms.trex_ctrl, "EN");
  auto en_aug = run_gvs_rho_sweep(
      rho_fn, preset.rho_grid_1d, cfg, ms.en_aug, ms.trex_ctrl, "EN+AUG");
  auto ien = run_gvs_rho_sweep(
      rho_fn, preset.rho_grid_1d, cfg, ms.ien, ms.trex_ctrl, "IEN");

  print_mc_rho_table(
      preset.scenario_tag,
      preset.rho_grid_1d,
      en,
      en_aug,
      ien,
      cfg,
      preset.file_stem_prefix + "_rho");
}

// =============================================================================
// Part 3 — 2-D SNR × rho sweep
// =============================================================================

/**
 * @brief Run a 2-D MC sweep over SNR and rho on the AR(1)-blocks DGP.
 */
static void run_ar1_part3_2d_sweep(const GVSSimConfig& cfg,
                                   const AR1BlocksPreset& preset,
                                   const AR1BlocksMethodSet& ms) {
  GVS2DDGPFactory dgp_2d = [&cfg](double snr, double rho, unsigned seed) {
    return make_ar1_blocks_dgp(cfg.n, cfg.p, snr, rho, seed);
  };

  cdiag::print_section_header("Part 3: 2-D SNR x rho | " + preset.scenario_tag);

  auto en = run_gvs_2d_sweep(
      dgp_2d, preset.snr_grid_2d, preset.rho_grid_2d, cfg,
      ms.en, ms.trex_ctrl, "EN");
  auto en_aug = run_gvs_2d_sweep(
      dgp_2d, preset.snr_grid_2d, preset.rho_grid_2d, cfg,
      ms.en_aug, ms.trex_ctrl, "EN+AUG");
  auto ien = run_gvs_2d_sweep(
      dgp_2d, preset.snr_grid_2d, preset.rho_grid_2d, cfg,
      ms.ien, ms.trex_ctrl, "IEN");

  std::vector<std::string> snr_labels, rho_labels;
  for (double s : preset.snr_grid_2d) {
    snr_labels.push_back("snr=" + std::to_string(s).substr(0, 4));
  }
  for (double r : preset.rho_grid_2d) {
    rho_labels.push_back("rho=" + std::to_string(r).substr(0, 4));
  }

  print_mc_matrix("mean_FDP [EN]", snr_labels, rho_labels, en, false);
  print_mc_matrix("mean_TPP [EN]", snr_labels, rho_labels, en, true);
  print_mc_matrix("mean_FDP [EN+AUG]", snr_labels, rho_labels, en_aug, false);
  print_mc_matrix("mean_TPP [EN+AUG]", snr_labels, rho_labels, en_aug, true);
  print_mc_matrix("mean_FDP [IEN]", snr_labels, rho_labels, ien, false);
  print_mc_matrix("mean_TPP [IEN]", snr_labels, rho_labels, ien, true);

  save_mc_2d_tables(
      preset.scenario_tag,
      snr_labels,
      rho_labels,
      en,
      en_aug,
      ien,
      cfg,
      preset.file_stem_prefix + "_2d");
}

// =============================================================================
// Top-level driver
// =============================================================================

/**
 * @brief Run the AR(1)-blocks benchmark through all three parts.
 */
static void run_ar1_blocks_benchmark(const GVSSimConfig& base_cfg,
                                     const AR1BlocksPreset& preset) {
  GVSSimConfig cfg = base_cfg;
  cfg.snr = preset.fixed_snr_for_rho;
  cfg.sd_x = 0.0; // not used; rho is passed directly to the AR(1) DGP

  const auto ms = make_method_set(cfg);

  cdiag::print_section_header(
      "AR(1)-within-block benchmark | " + preset.scenario_tag +
      "\nDGP = make_ar1_blocks_dgp"
      "\nblocks = {20, 50, 80, 65}, active = {0,1,2}, inactive trap = {3}");

  run_ar1_part1_snr_sweep(cfg, preset, ms);
  run_ar1_part2_rho_sweep(cfg, preset, ms);
  run_ar1_part3_2d_sweep(cfg, preset, ms);

  std::cout << preset.scenario_tag << " GVS MC simulations complete.\n";
}

// =============================================================================
// main
// =============================================================================

int main() {

  std::cout.setf(std::ios::unitbuf);
  omp_set_num_threads(6);
  std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

  GVSSimConfig cfg;
  cfg.n = 200;
  cfg.p = 500;
  cfg.sd_x = 0.0; // unused in this file
  cfg.tFDR = 0.1;
  cfg.K = 20;
  cfg.num_MC = 200;
  cfg.base_seed = 2026;
  cfg.corr_max = 0.98;
  cfg.snr = 2.0;

  run_ar1_blocks_benchmark(cfg, make_ar1_blocks_preset());

  std::cout << "All AR(1)-blocks benchmark simulations complete.\n";
  return 0;
}
