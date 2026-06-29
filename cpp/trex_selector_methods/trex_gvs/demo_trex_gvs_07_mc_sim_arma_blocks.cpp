// ==============================================================================
// demo_trex_gvs_07_mc_sim_arma_blocks.cpp
// ==============================================================================
/**
 * @file demo_trex_gvs_07_mc_sim_arma_blocks.cpp
 *
 * @brief MC simulation: T-Rex+GVS on heterogeneous ARMA blocks.
 *
 * @details
 * Three parts:
 *   1. SNR sweep at fixed ar_coef = 0.8
 *   2. ar_coef sweep at fixed SNR = 2.0
 *   3. 2-D SNR × ar_coef sweep
 *
 * DGP:
 *   make_arma_blocks_dgp()
 *   - 4 blocks of sizes {20, 50, 80, 65}
 *   - heterogeneous within-block ARMA structure
 *   - blocks 0, 1, 2 active with beta = 3
 *   - block 3 is an inactive AR(1) trap
 *   - support size s = 150
 *
 * Shared selector settings:
 *   - K = 20
 *   - tFDR = 0.1
 *   - corr_max = 0.98
 *   - single-linkage HAC
 *   - lambda2 = CV_1SE_CCD
 *
 * Methods compared:
 *   - EN (TENET)
 *   - EN (TENET_AUG)
 *   - IEN
 *
 * MC:
 *   - 200 trials per grid point
 */
// ==============================================================================

#include "trex_gvs_dgps.hpp"
#include "trex_gvs_simulation_utils.hpp"

// ==============================================================================

using namespace gvs_demo;

// =============================================================================
// Preset and method bundle
// =============================================================================

struct ARMABlocksPreset {
    std::string scenario_tag;
    std::string file_stem_prefix;

    double fixed_ar_coef_for_snr = 0.8;
    double fixed_snr_for_ar_coef = 2.0;

    std::vector<double> snr_grid_1d;
    std::vector<double> ar_coef_grid_1d;
    std::vector<double> snr_grid_2d;
    std::vector<double> ar_coef_grid_2d;
};

struct ARMABlocksMethodSet {
    TRexControlParameter trex_ctrl;
    TRexGVSControlParameter en;
    TRexGVSControlParameter en_aug;
    TRexGVSControlParameter ien;
};

// =============================================================================
// Helpers
// =============================================================================

static ARMABlocksPreset make_arma_blocks_preset()
{
    ARMABlocksPreset p;
    p.scenario_tag = "ARMA-Blocks";
    p.file_stem_prefix = "gvs_arma_blocks";

    p.fixed_ar_coef_for_snr = 0.8;
    p.fixed_snr_for_ar_coef = 2.0;

    p.snr_grid_1d = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    p.ar_coef_grid_1d = {
        0.10, 0.20, 0.30, 0.40,
        0.50, 0.60, 0.70, 0.80,
        0.90, 0.95
    };

    p.snr_grid_2d = {0.2, 0.5, 1.0, 2.0, 5.0};

    p.ar_coef_grid_2d = {
        0.30, 0.50, 0.70, 0.80,
        0.90, 0.95
    };

    return p;
}

static ARMABlocksMethodSet make_method_set(const GVSSimConfig& cfg)
{
    ARMABlocksMethodSet ms;
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
 * @brief Run an MC SNR sweep on the ARMA-blocks DGP at fixed ar_coef.
 */
static void run_arma_part1_snr_sweep(
    const GVSSimConfig& cfg,
    const ARMABlocksPreset& preset,
    const ARMABlocksMethodSet& ms)
{
    const double fixed_ar = preset.fixed_ar_coef_for_snr;

    GVSDGPFactory snr_fn = [&cfg, fixed_ar](double snr, unsigned seed) {
        return make_arma_blocks_dgp(cfg.n, cfg.p, snr, fixed_ar, seed);
    };

    cdiag::print_section_header(
        "Part 1: SNR Sweep | " + preset.scenario_tag +
        "\nn=" + std::to_string(cfg.n) +
        " p=" + std::to_string(cfg.p) +
        " s=150 MC=" + std::to_string(cfg.num_MC) +
        " tFDR=" + std::to_string(cfg.tFDR).substr(0, 3) +
        " ar_coef=0.8 (fixed)");

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
// Part 2 — ar_coef sweep
// =============================================================================

/**
 * @brief Run an MC ar_coef sweep on the ARMA-blocks DGP at fixed SNR.
 */
static void run_arma_part2_arcoef_sweep(
    const GVSSimConfig& cfg,
    const ARMABlocksPreset& preset,
    const ARMABlocksMethodSet& ms)
{
    // GVSRhoDGPFactory is reused here as a generic 1-D sweep factory:
    // the first argument is ar_coef, not rho.
    GVSRhoDGPFactory ar_fn = [&cfg](double ar_coef, unsigned seed) {
        return make_arma_blocks_dgp(cfg.n, cfg.p, cfg.snr, ar_coef, seed);
    };

    cdiag::print_section_header(
        "Part 2: ar_coef Sweep | " + preset.scenario_tag +
        " (SNR=" + std::to_string(cfg.snr).substr(0, 3) + ")");

    auto en = run_gvs_rho_sweep(
        ar_fn, preset.ar_coef_grid_1d, cfg, ms.en, ms.trex_ctrl, "EN");

    auto en_aug = run_gvs_rho_sweep(
        ar_fn, preset.ar_coef_grid_1d, cfg, ms.en_aug, ms.trex_ctrl, "EN+AUG");

    auto ien = run_gvs_rho_sweep(
        ar_fn, preset.ar_coef_grid_1d, cfg, ms.ien, ms.trex_ctrl, "IEN");

    print_mc_param_sweep_table(
        preset.scenario_tag,
        "ar_coef",
        preset.ar_coef_grid_1d,
        en,
        en_aug,
        ien,
        cfg,
        preset.file_stem_prefix + "_arcoef");
}

// =============================================================================
// Part 3 — 2-D SNR × ar_coef sweep
// =============================================================================

/**
 * @brief Run a 2-D MC sweep over SNR and ar_coef on the ARMA-blocks DGP.
 */
static void run_arma_part3_2d_sweep(
    const GVSSimConfig& cfg,
    const ARMABlocksPreset& preset,
    const ARMABlocksMethodSet& ms)
{
    GVS2DDGPFactory dgp_2d = [&cfg](double snr, double ar_coef, unsigned seed) {
        return make_arma_blocks_dgp(cfg.n, cfg.p, snr, ar_coef, seed);
    };

    cdiag::print_section_header(
        "Part 3: 2-D SNR x ar_coef | " + preset.scenario_tag);

    auto en = run_gvs_2d_sweep(
        dgp_2d,
        preset.snr_grid_2d,
        preset.ar_coef_grid_2d,
        cfg,
        ms.en,
        ms.trex_ctrl,
        "EN");

    auto en_aug = run_gvs_2d_sweep(
        dgp_2d,
        preset.snr_grid_2d,
        preset.ar_coef_grid_2d,
        cfg,
        ms.en_aug,
        ms.trex_ctrl,
        "EN+AUG");

    auto ien = run_gvs_2d_sweep(
        dgp_2d,
        preset.snr_grid_2d,
        preset.ar_coef_grid_2d,
        cfg,
        ms.ien,
        ms.trex_ctrl,
        "IEN");

    std::vector<std::string> snr_labels, ar_labels;
    for (double s : preset.snr_grid_2d)
        snr_labels.push_back("snr=" + std::to_string(s).substr(0, 4));
    for (double a : preset.ar_coef_grid_2d)
        ar_labels.push_back("ac=" + std::to_string(a).substr(0, 4));

    print_mc_matrix("mean_FDP [EN]",     snr_labels, ar_labels, en,     false);
    print_mc_matrix("mean_TPP [EN]",     snr_labels, ar_labels, en,     true);
    print_mc_matrix("mean_FDP [EN+AUG]", snr_labels, ar_labels, en_aug, false);
    print_mc_matrix("mean_TPP [EN+AUG]", snr_labels, ar_labels, en_aug, true);
    print_mc_matrix("mean_FDP [IEN]",    snr_labels, ar_labels, ien,    false);
    print_mc_matrix("mean_TPP [IEN]",    snr_labels, ar_labels, ien,    true);

    save_mc_2d_tables(
        preset.scenario_tag,
        snr_labels,
        ar_labels,
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
 * @brief Run the heterogeneous ARMA-blocks benchmark through all three parts.
 */
static void run_arma_blocks_benchmark(
    const GVSSimConfig& base_cfg,
    const ARMABlocksPreset& preset)
{
    GVSSimConfig cfg = base_cfg;
    cfg.snr = preset.fixed_snr_for_ar_coef;

    const auto ms = make_method_set(cfg);

    cdiag::print_section_header(
        "ARMA-blocks benchmark | " + preset.scenario_tag +
        "\nDGP = make_arma_blocks_dgp"
        "\nblocks = {20, 50, 80, 65}, s = 150"
        "\nactive blocks = {0, 1, 2}, trap block = 3");

    run_arma_part1_snr_sweep(cfg, preset, ms);
    run_arma_part2_arcoef_sweep(cfg, preset, ms);
    run_arma_part3_2d_sweep(cfg, preset, ms);

    std::cout << preset.scenario_tag << " GVS MC simulations complete.\n";
}

// =============================================================================
// main
// =============================================================================

int main()
{
    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    GVSSimConfig cfg;
    cfg.n         = 200;
    cfg.p         = 500;
    cfg.sd_x      = 0.0;   // not used in this DGP; ar_coef is passed directly
    cfg.tFDR      = 0.1;
    cfg.K         = 20;
    cfg.num_MC    = 200;
    cfg.base_seed = 2026;
    cfg.corr_max  = 0.98;
    cfg.snr       = 2.0;   // fixed Part-2 value; overwritten in driver

    run_arma_blocks_benchmark(cfg, make_arma_blocks_preset());

    std::cout << "All ARMA-block benchmark simulations complete.\n";
    return 0;
}
