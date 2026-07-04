// ==============================================================================
// demo_trex_gvs_02_mc_sim_scattered_blocks.cpp
// ==============================================================================
/**
 * @file demo_trex_gvs_02_mc_sim_scattered_blocks.cpp
 *
 * @brief MC simulation: T-Rex+GVS on the scattered-grouped DGP.
 *
 * @details
 * This scenario evaluates grouped-signal recovery when active variables are
 * randomly scattered across the design rather than arranged in contiguous
 * blocks. Within-group correlation is controlled through sd_x via
 * rho = 1 / (1 + sd_x^2).
 *
 * Three parts:
 * 1. SNR sweep at fixed sd_x = 0.1 (equivalently sqrt(0.01), rho ≈ 0.99)
 * 2. rho sweep at fixed SNR = 2.0, using sd_x = sqrt((1 - rho)/rho)
 * 3. 2-D SNR × rho sweep
 *
 * DGP:
 * - make_scattered_grouped_dgp()
 * - n_groups = 3
 * - group_size = 50
 * - all grouped variables active
 * - support size s = 150
 * - p = 500 total variables, with remaining columns i.i.d. noise
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

struct ScatteredGroupedPreset {
  std::string scenario_tag;
  std::string file_stem_prefix;

  int n_groups = 3;
  int group_size = 50;

  double fixed_sd_x_for_snr = 0.1; // sqrt(0.01)
  double fixed_snr_for_rho = 2.0;

  std::vector<double> snr_grid_1d;
  std::vector<double> rho_grid_1d;
  std::vector<double> snr_grid_2d;
  std::vector<double> rho_grid_2d;
};

struct ScatteredGroupedMethodSet {
  TRexControlParameter trex_ctrl;
  TRexGVSControlParameter en;
  TRexGVSControlParameter en_aug;
  TRexGVSControlParameter ien;
};

// =============================================================================
// Helpers
// =============================================================================

static double rho_to_sd_x(double rho) {
  return std::sqrt((1.0 - rho) / rho);
}

static ScatteredGroupedPreset make_scattered_grouped_preset() {
  ScatteredGroupedPreset p;
  p.scenario_tag = "Scattered-Grouped";
  p.file_stem_prefix = "gvs_scattered_grouped";

  p.n_groups = 3;
  p.group_size = 50;

  p.fixed_sd_x_for_snr = 0.1;
  p.fixed_snr_for_rho = 2.0;

  p.snr_grid_1d = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_1d = {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99};
  p.snr_grid_2d = {0.2, 0.5, 1.0, 2.0, 5.0};
  p.rho_grid_2d = {0.30, 0.50, 0.70, 0.90, 0.95, 0.99};

  return p;
}

static ScatteredGroupedMethodSet make_method_set(const GVSSimConfig& cfg) {
  ScatteredGroupedMethodSet ms;
  ms.trex_ctrl = make_gvs_trex_control(cfg.K);

  ms.en.gvs_type = GVSType::EN;
  ms.en.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
  ms.en.corr_max = cfg.corr_max;
  ms.en.hc_linkage = hac::LinkageMethod::Single;
  ms.en.en_solver = ENSolverType::TENET;

  ms.en_aug = ms.en; // copy and overwrite solver type
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
 * @brief Run an MC SNR sweep on the scattered-grouped DGP at fixed sd_x.
 */
static void run_scattered_part1_snr_sweep(
    const GVSSimConfig& cfg,
    const ScatteredGroupedPreset& preset,
    const ScatteredGroupedMethodSet& ms) {
  const int n_groups = preset.n_groups;
  const int group_size = preset.group_size;
  const double fixed_sd_x = preset.fixed_sd_x_for_snr;

  GVSDGPFactory snr_fn = [&cfg, n_groups, group_size, fixed_sd_x](
                             double snr, unsigned seed) {
    return make_scattered_grouped_dgp(
        cfg.n, cfg.p, n_groups, group_size, snr, fixed_sd_x, seed);
  };

  cdiag::print_section_header(
      "Part 1: SNR Sweep | " + preset.scenario_tag +
      "\nn=" + std::to_string(cfg.n) + " p=" + std::to_string(cfg.p) +
      " n_groups=" + std::to_string(n_groups) +
      " group_size=" + std::to_string(group_size) +
      " s=150 MC=" + std::to_string(cfg.num_MC) +
      " tFDR=" + fmt_num(cfg.tFDR) +
      " sd_x=0.1 (fixed)");

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
 * @brief Run an MC rho sweep on the scattered-grouped DGP at fixed SNR.
 */
static void run_scattered_part2_rho_sweep(
    const GVSSimConfig& cfg,
    const ScatteredGroupedPreset& preset,
    const ScatteredGroupedMethodSet& ms) {
  const int n_groups = preset.n_groups;
  const int group_size = preset.group_size;

  GVSRhoDGPFactory rho_fn = [&cfg, n_groups, group_size](
                                double rho, unsigned seed) {
    const double sd_x = rho_to_sd_x(rho);
    return make_scattered_grouped_dgp(
        cfg.n, cfg.p, n_groups, group_size, cfg.snr, sd_x, seed);
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
 * @brief Run a 2-D MC sweep over SNR and rho on the scattered-grouped DGP.
 */
static void run_scattered_part3_2d_sweep(
    const GVSSimConfig& cfg,
    const ScatteredGroupedPreset& preset,
    const ScatteredGroupedMethodSet& ms) {
  const int n_groups = preset.n_groups;
  const int group_size = preset.group_size;

  GVS2DDGPFactory dgp_2d = [&cfg, n_groups, group_size](
                               double snr, double rho, unsigned seed) {
    const double sd_x = rho_to_sd_x(rho);
    return make_scattered_grouped_dgp(
        cfg.n, cfg.p, n_groups, group_size, snr, sd_x, seed);
  };

  cdiag::print_section_header(
      "Part 3: 2-D SNR x rho | " + preset.scenario_tag);

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

  print_mc_matrix("mean_TPP [EN]", snr_labels, rho_labels, en, true);
  print_mc_matrix("mean_FDP [EN]", snr_labels, rho_labels, en, false);
  print_mc_matrix("mean_TPP [EN+AUG]", snr_labels, rho_labels, en_aug, true);
  print_mc_matrix("mean_FDP [EN+AUG]", snr_labels, rho_labels, en_aug, false);
  print_mc_matrix("mean_TPP [IEN]", snr_labels, rho_labels, ien, true);
  print_mc_matrix("mean_FDP [IEN]", snr_labels, rho_labels, ien, false);

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
 * @brief Run the scattered-grouped benchmark through all three parts.
 */
static void run_scattered_grouped_benchmark(
    const GVSSimConfig& base_cfg,
    const ScatteredGroupedPreset& preset) {
  GVSSimConfig cfg = base_cfg;
  cfg.sd_x = preset.fixed_sd_x_for_snr;
  cfg.snr = preset.fixed_snr_for_rho;

  const auto ms = make_method_set(cfg);

  cdiag::print_section_header(
      "Scattered-grouped benchmark | " + preset.scenario_tag +
      "\nDGP = make_scattered_grouped_dgp"
      "\nn_groups = " + std::to_string(preset.n_groups) +
      ", group_size = " + std::to_string(preset.group_size) +
      ", s = " + std::to_string(preset.n_groups * preset.group_size));

  run_scattered_part1_snr_sweep(cfg, preset, ms);
  run_scattered_part2_rho_sweep(cfg, preset, ms);
  run_scattered_part3_2d_sweep(cfg, preset, ms);

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
  cfg.sd_x = 0.1;   // fixed Part-1 value; overwritten again in top-level driver
  cfg.tFDR = 0.1;
  cfg.K = 20;
  cfg.num_MC = 200;
  cfg.base_seed = 2026;
  cfg.corr_max = 0.98;
  cfg.snr = 2.0;    // fixed Part-2 value; overwritten again in top-level driver


  // ==========================================================================
  //  Run simulations
  // ==========================================================================

  if (true) run_scattered_grouped_benchmark(cfg, make_scattered_grouped_preset());

  std::cout << "All scattered-grouped benchmark simulations complete.\n";
  return 0;
}
