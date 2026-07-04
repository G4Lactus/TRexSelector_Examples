// ==============================================================================
// demo_trex_gvs_03_mc_sim_mixed_blocks.cpp
// ==============================================================================
/**
 * @file demo_trex_gvs_03_mc_sim_mixed_blocks.cpp
 *
 * @brief MC simulation: T-Rex+GVS on the mixed-blocks DGP, with two benchmark
 * presets sharing the same underlying generator.
 *
 * @details
 * This file merges the former "mixed-blocks" and "equicorr-blocks" drivers into
 * a single implementation. Both presets use make_mixed_blocks_dgp(), i.e. a
 * four-block Gaussian design with block sizes {20, 50, 80, 65}, three active
 * blocks (total support size s = 150), and one inactive equicorrelated trap
 * block, placed in random order with random noise gaps.
 *
 * Presets:
 * 1. Mixed-Blocks
 *    - Part 1: SNR sweep at fixed sd_x = sqrt(0.01)   (rho ≈ 0.99)
 *    - Part 2: rho sweep at fixed SNR = 2.0
 *    - Part 3: 2-D SNR × rho grid
 *
 * 2. Mixed-Blocks-Rho075
 *    - Part 1: SNR sweep at fixed rho = 0.75
 *    - Part 2: rho sweep at fixed SNR = 2.0
 *    - Part 3: 2-D SNR × rho grid with slightly expanded rho coverage
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
// Benchmark presets
// =============================================================================

enum class MixedBlocksMode {
  DEFAULT_SD_X,
  FIXED_RHO_075
};

struct MixedBlocksPreset {
  MixedBlocksMode mode;
  std::string scenario_tag;
  std::string file_stem_prefix;

  double fixed_sd_x_for_snr = 0.1;
  double fixed_rho_for_snr = 0.75;
  double fixed_snr_for_rho = 2.0;

  std::vector<double> snr_grid_1d;
  std::vector<double> rho_grid_1d;
  std::vector<double> snr_grid_2d;
  std::vector<double> rho_grid_2d;
};

struct MixedBlocksMethodSet {
  TRexControlParameter trex_ctrl;
  TRexGVSControlParameter en;
  TRexGVSControlParameter en_aug;
  TRexGVSControlParameter ien;
};

// =============================================================================
// Preset builders
// =============================================================================

static MixedBlocksPreset make_default_sd_x_preset() {
  MixedBlocksPreset p;
  p.mode = MixedBlocksMode::DEFAULT_SD_X;
  p.scenario_tag = "Mixed-Blocks";
  p.file_stem_prefix = "gvs_mixed_blocks";
  p.fixed_sd_x_for_snr = std::sqrt(0.01);
  p.fixed_snr_for_rho = 2.0;
  p.snr_grid_1d = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_1d = {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99};
  p.snr_grid_2d = {0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_2d = {0.30, 0.50, 0.70, 0.90, 0.95, 0.99};
  return p;
}

static MixedBlocksPreset make_fixed_rho_preset() {
  MixedBlocksPreset p;
  p.mode = MixedBlocksMode::FIXED_RHO_075;
  p.scenario_tag = "Mixed-Blocks-Rho075";
  p.file_stem_prefix = "gvs_mixed_blocks_rho075";
  p.fixed_rho_for_snr = 0.75;
  p.fixed_snr_for_rho = 2.0;
  p.snr_grid_1d = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_1d = {0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99};
  p.snr_grid_2d = {0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_2d = {0.30, 0.50, 0.70, 0.85, 0.90, 0.99};
  return p;
}

// =============================================================================
// Shared method controls
// =============================================================================

static MixedBlocksMethodSet make_method_set(const GVSSimConfig& cfg) {
  MixedBlocksMethodSet ms;
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
// Helpers
// =============================================================================

static double resolve_fixed_sd_x_for_snr(const MixedBlocksPreset& preset) {
  if (preset.mode == MixedBlocksMode::DEFAULT_SD_X) {
    return preset.fixed_sd_x_for_snr;
  }
  return std::sqrt((1.0 - preset.fixed_rho_for_snr) / preset.fixed_rho_for_snr);
}

static std::string make_part1_header_suffix(const MixedBlocksPreset& preset) {
  if (preset.mode == MixedBlocksMode::DEFAULT_SD_X) {
    return "sd_x=sqrt(0.01)";
  }
  return "rho=0.75 (fixed)";
}

// =============================================================================
// Part 1 — SNR sweep
// =============================================================================

/**
 * @brief Run an MC SNR sweep on the mixed-blocks DGP using the selected preset.
 *
 * For DEFAULT_SD_X, fixes sd_x = sqrt(0.01). For FIXED_RHO_075, fixes rho=0.75
 * and converts to sd_x = sqrt((1-rho)/rho). Compares EN (TENET), EN (TENET_AUG),
 * and IEN across the SNR grid.
 */
static void run_part1_snr_sweep(const GVSSimConfig& cfg,
                                const MixedBlocksPreset& preset,
                                const MixedBlocksMethodSet& ms) {
  const double fixed_sd_x = resolve_fixed_sd_x_for_snr(preset);

  GVSDGPFactory snr_fn = [&cfg, fixed_sd_x](double snr, unsigned seed) {
    return make_mixed_blocks_dgp(cfg.n, cfg.p, snr, fixed_sd_x, seed);
  };

  cdiag::print_section_header(
      "Part 1: SNR Sweep | " + preset.scenario_tag +
      "\nn=" + std::to_string(cfg.n) + " p=" + std::to_string(cfg.p) +
      " s=150 MC=" + std::to_string(cfg.num_MC) +
      " tFDR=" + fmt_num(cfg.tFDR) + " " +
      make_part1_header_suffix(preset));

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
 * @brief Run an MC rho sweep on the mixed-blocks DGP at fixed SNR.
 *
 * Converts each target rho to sd_x = sqrt((1 - rho) / rho), then compares
 * EN (TENET), EN (TENET_AUG), and IEN across the rho grid.
 */
static void run_part2_rho_sweep(const GVSSimConfig& cfg,
                                const MixedBlocksPreset& preset,
                                const MixedBlocksMethodSet& ms) {
  GVSRhoDGPFactory rho_fn = [&cfg](double rho, unsigned seed) {
    const double sd_x = std::sqrt((1.0 - rho) / rho);
    return make_mixed_blocks_dgp(cfg.n, cfg.p, cfg.snr, sd_x, seed);
  };

  cdiag::print_section_header(
      "Part 2: rho Sweep | " + preset.scenario_tag +
      " (SNR=" + fmt_num(cfg.snr) + ")");

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
 * @brief Run a 2-D MC sweep over SNR and rho on the mixed-blocks DGP.
 *
 * Evaluates EN (TENET), EN (TENET_AUG), and IEN on an SNR × rho grid, then
 * prints and saves matrix summaries for mean FDP and mean TPP.
 */
static void run_part3_2d_sweep(const GVSSimConfig& cfg,
                               const MixedBlocksPreset& preset,
                               const MixedBlocksMethodSet& ms) {
  GVS2DDGPFactory dgp_2d = [&cfg](double snr, double rho, unsigned seed) {
    const double sd_x = std::sqrt((1.0 - rho) / rho);
    return make_mixed_blocks_dgp(cfg.n, cfg.p, snr, sd_x, seed);
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
    snr_labels.push_back("snr=" + fmt_num(s));
  }
  for (double r : preset.rho_grid_2d) {
    rho_labels.push_back("rho=" + fmt_num(r));
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
// Top-level benchmark driver
// =============================================================================

/**
 * @brief Run one mixed-blocks benchmark preset through all three parts.
 */
static void run_mixed_blocks_benchmark(const GVSSimConfig& base_cfg,
                                       const MixedBlocksPreset& preset) {
  GVSSimConfig cfg = base_cfg;
  cfg.snr = preset.fixed_snr_for_rho;

  const auto ms = make_method_set(cfg);

  cdiag::print_section_header(
      "Mixed-blocks benchmark | " + preset.scenario_tag +
      "\nDGP = make_mixed_blocks_dgp"
      "\nblocks = {20, 50, 80, 65}, active = {0,1,2}, inactive trap = {3}");

  run_part1_snr_sweep(cfg, preset, ms);
  run_part2_rho_sweep(cfg, preset, ms);
  run_part3_2d_sweep(cfg, preset, ms);

  std::cout << preset.scenario_tag << " GVS MC simulations complete.\n";
}

// =============================================================================
// main
// =============================================================================

int main() {

  // ==========================================================================
  //  Simulation parameters
  // ==========================================================================

  std::cout.setf(std::ios::unitbuf);
  omp_set_num_threads(6);
  std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

  GVSSimConfig cfg;
  cfg.n = 200;
  cfg.p = 500;
  cfg.tFDR = 0.1;
  cfg.K = 20;
  cfg.num_MC = 200;
  cfg.base_seed = 2026;
  cfg.corr_max = 0.98;
  cfg.snr = 2.0;

  // ==========================================================================
  //  Run simulations
  // ==========================================================================

  if (true) run_mixed_blocks_benchmark(cfg, make_default_sd_x_preset());
  if (true) run_mixed_blocks_benchmark(cfg, make_fixed_rho_preset());

  std::cout << "All mixed-blocks benchmark presets complete.\n";
  return 0;
}
