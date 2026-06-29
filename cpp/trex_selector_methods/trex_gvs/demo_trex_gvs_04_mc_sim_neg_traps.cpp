// ==============================================================================
// demo_trex_gvs_04_mc_sim_neg_traps.cpp
// ==============================================================================
/**
 * @file demo_trex_gvs_04_mc_sim_neg_traps.cpp
 *
 * @brief MC simulation: T-Rex+GVS on the negative-traps DGP.
 *  This scenario evaluates recovery of an active negatively correlated group
 *  in the presence of two spatially separated inactive trap groups with the
 *  same sign-flipped correlation structure.
 *
 * @details
 *  Three parts, each in its own function:
 *    sim_neg_traps_part1() — SNR sweep (fixed sd_x = sqrt(0.01))
 *    sim_neg_traps_part2() — rho sweep (fixed SNR = 2.0, sd_x derived from rho)
 *    sim_neg_traps_part3() — 2-D SNR × rho grid
 *
 *  DGP: active group (0–99) + Trap1 (100–199) + noise (200–299) +
 *       Trap2 (300–359) + noise (360+); s = 100.
 *  In each correlated group, one half loads positively and the other half
 *  negatively on a latent factor, inducing strong within-group sign structure.
 *  Selector: K=20, tFDR=0.1, corr_max=0.98, lambda2=CV_1SE_CCD.
 *  Methods compared: EN (TENET), EN (TENET_AUG), IEN.
 *  MC: 200 trials per grid point.
 */
// ==============================================================================

// project includes
#include "trex_gvs_dgps.hpp"
#include "trex_gvs_simulation_utils.hpp"

// ==============================================================================

using namespace gvs_demo;

// ==============================================================================
//  Part 1 — SNR sweep
// ==============================================================================

/**
 * @brief Run an MC SNR sweep on the negative-traps DGP at fixed within-group correlation.
 *
 * Uses sd_x = sqrt(0.01) (rho ≈ 0.99) and compares EN (TENET), EN (TENET_AUG),
 * and IEN across the SNR grid; prints and saves mean FDP/TPP tables.
 */

static void sim_neg_traps_part1(const GVSSimConfig& cfg, const std::string& tag)
{
    const std::vector<double> snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    GVSDGPFactory snr_fn = [&cfg](double snr, unsigned seed) {
        return make_neg_traps_dgp(cfg.n, cfg.p, snr, cfg.sd_x, seed);
    };

    auto trex_ctrl = make_gvs_trex_control(cfg.K);

    TRexGVSControlParameter gvs_ctrl_en;
    gvs_ctrl_en.gvs_type       = GVSType::EN;
    gvs_ctrl_en.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_en.corr_max       = cfg.corr_max;
    gvs_ctrl_en.hc_linkage     = hac::LinkageMethod::Single;
    gvs_ctrl_en.en_solver      = ENSolverType::TENET;

    TRexGVSControlParameter gvs_ctrl_en_aug;
    gvs_ctrl_en_aug.gvs_type       = GVSType::EN;
    gvs_ctrl_en_aug.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_en_aug.corr_max       = cfg.corr_max;
    gvs_ctrl_en_aug.hc_linkage     = hac::LinkageMethod::Single;
    gvs_ctrl_en_aug.en_solver      = ENSolverType::TENET_AUG;

    TRexGVSControlParameter gvs_ctrl_ien;
    gvs_ctrl_ien.gvs_type       = GVSType::IEN;
    gvs_ctrl_ien.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_ien.corr_max       = cfg.corr_max;
    gvs_ctrl_ien.hc_linkage     = hac::LinkageMethod::Single;

    cdiag::print_section_header(
        "Part 1: SNR Sweep  |  " + tag +
        "\nn=" + std::to_string(cfg.n) + "  p=" + std::to_string(cfg.p) +
        "  s=100  MC=" + std::to_string(cfg.num_MC) +
        "  tFDR=" + std::to_string(cfg.tFDR).substr(0, 3) +
        "  sd_x=sqrt(0.01)");

    auto en_snr     = run_gvs_snr_sweep(
        snr_fn, snr_grid, cfg, gvs_ctrl_en,     trex_ctrl, "EN");

    auto en_aug_snr = run_gvs_snr_sweep(
        snr_fn, snr_grid, cfg, gvs_ctrl_en_aug, trex_ctrl, "EN+AUG");

    auto ien_snr    = run_gvs_snr_sweep(
        snr_fn, snr_grid, cfg, gvs_ctrl_ien,    trex_ctrl, "IEN");

    print_mc_snr_table(tag, snr_grid, en_snr, en_aug_snr, ien_snr, cfg, "gvs_" + tag + "_snr");
}


// ==============================================================================
//  Part 2 — rho sweep
// ==============================================================================

/**
 * @brief Run an MC rho sweep on the negative-traps DGP at fixed SNR.
 *
 * Converts each target rho to sd_x = sqrt((1 - rho) / rho), then compares
 * EN (TENET), EN (TENET_AUG), and IEN across the rho grid.
 */

static void sim_neg_traps_part2(const GVSSimConfig& cfg, const std::string& tag)
{
    const std::vector<double> rho_grid = {0.10, 0.20, 0.30, 0.50, 0.70,
                                           0.80, 0.90, 0.95, 0.99};

    GVSRhoDGPFactory rho_fn = [&cfg](double rho, unsigned seed) {
        const double sd_x = std::sqrt((1.0 - rho) / rho);
        return make_neg_traps_dgp(cfg.n, cfg.p, cfg.snr, sd_x, seed);
    };

    auto trex_ctrl = make_gvs_trex_control(cfg.K);

    TRexGVSControlParameter gvs_ctrl_en;
    gvs_ctrl_en.gvs_type       = GVSType::EN;
    gvs_ctrl_en.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_en.corr_max       = cfg.corr_max;
    gvs_ctrl_en.hc_linkage     = hac::LinkageMethod::Single;
    gvs_ctrl_en.en_solver      = ENSolverType::TENET;

    TRexGVSControlParameter gvs_ctrl_en_aug;
    gvs_ctrl_en_aug.gvs_type       = GVSType::EN;
    gvs_ctrl_en_aug.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_en_aug.corr_max       = cfg.corr_max;
    gvs_ctrl_en_aug.hc_linkage     = hac::LinkageMethod::Single;
    gvs_ctrl_en_aug.en_solver      = ENSolverType::TENET_AUG;

    TRexGVSControlParameter gvs_ctrl_ien;
    gvs_ctrl_ien.gvs_type       = GVSType::IEN;
    gvs_ctrl_ien.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_ien.corr_max       = cfg.corr_max;
    gvs_ctrl_ien.hc_linkage     = hac::LinkageMethod::Single;

    cdiag::print_section_header(
        "Part 2: rho Sweep  |  " + tag +
        "  (SNR=" + std::to_string(cfg.snr).substr(0, 3) + ")");

    auto en_rho     = run_gvs_rho_sweep(
        rho_fn, rho_grid, cfg, gvs_ctrl_en,     trex_ctrl, "EN");
    auto en_aug_rho = run_gvs_rho_sweep(
        rho_fn, rho_grid, cfg, gvs_ctrl_en_aug, trex_ctrl, "EN+AUG");
    auto ien_rho    = run_gvs_rho_sweep(
        rho_fn, rho_grid, cfg, gvs_ctrl_ien,    trex_ctrl, "IEN");
    print_mc_rho_table(tag, rho_grid, en_rho, en_aug_rho, ien_rho, cfg, "gvs_" + tag + "_rho");
}


// ==============================================================================
//  Part 3 — 2-D SNR × rho sweep
// ==============================================================================

/**
 * @brief Run a 2-D MC sweep over SNR and rho on the negative-traps DGP.
 *
 * Evaluates EN (TENET), EN (TENET_AUG), and IEN on an SNR × rho grid, then
 * prints and saves matrix summaries for mean TPP and mean FDP.
 */

static void sim_neg_traps_part3(const GVSSimConfig& cfg, const std::string& tag)
{
    const std::vector<double> snr_grid_2d = {0.2, 0.5, 1.0, 2.0, 5.0};
    const std::vector<double> rho_grid_2d = {
        0.30, 0.50, 0.70, 0.90, 0.95, 0.99};

    GVS2DDGPFactory dgp_2d = [&cfg](double snr, double rho, unsigned seed) {
        const double sd_x = std::sqrt((1.0 - rho) / rho);
        return make_neg_traps_dgp(cfg.n, cfg.p, snr, sd_x, seed);
    };

    auto trex_ctrl = make_gvs_trex_control(cfg.K);

    TRexGVSControlParameter gvs_ctrl_en;
    gvs_ctrl_en.gvs_type       = GVSType::EN;
    gvs_ctrl_en.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_en.corr_max       = cfg.corr_max;
    gvs_ctrl_en.hc_linkage     = hac::LinkageMethod::Single;
    gvs_ctrl_en.en_solver      = ENSolverType::TENET;

    TRexGVSControlParameter gvs_ctrl_en_aug;
    gvs_ctrl_en_aug.gvs_type       = GVSType::EN;
    gvs_ctrl_en_aug.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_en_aug.corr_max       = cfg.corr_max;
    gvs_ctrl_en_aug.hc_linkage     = hac::LinkageMethod::Single;
    gvs_ctrl_en_aug.en_solver      = ENSolverType::TENET_AUG;

    TRexGVSControlParameter gvs_ctrl_ien;
    gvs_ctrl_ien.gvs_type       = GVSType::IEN;
    gvs_ctrl_ien.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
    gvs_ctrl_ien.corr_max       = cfg.corr_max;
    gvs_ctrl_ien.hc_linkage     = hac::LinkageMethod::Single;

    cdiag::print_section_header("Part 3: 2-D SNR x rho  |  " + tag);

    auto en_2d     = run_gvs_2d_sweep(dgp_2d, snr_grid_2d, rho_grid_2d, cfg,
                                       gvs_ctrl_en,     trex_ctrl, "EN");
    auto en_aug_2d = run_gvs_2d_sweep(dgp_2d, snr_grid_2d, rho_grid_2d, cfg,
                                       gvs_ctrl_en_aug, trex_ctrl, "EN+AUG");
    auto ien_2d    = run_gvs_2d_sweep(dgp_2d, snr_grid_2d, rho_grid_2d, cfg,
                                       gvs_ctrl_ien,    trex_ctrl, "IEN");

    std::vector<std::string> snr_labels, rho_labels;
    for (double s : snr_grid_2d)
        snr_labels.push_back("snr=" + std::to_string(s).substr(0, 4));
    for (double r : rho_grid_2d)
        rho_labels.push_back("rho=" + std::to_string(r).substr(0, 4));

    print_mc_matrix("mean_TPP [EN]",     snr_labels,
                    rho_labels, en_2d,     true);
    print_mc_matrix("mean_FDP [EN]",     snr_labels,
                    rho_labels, en_2d,     false);
    print_mc_matrix("mean_TPP [EN+AUG]", snr_labels,
                    rho_labels, en_aug_2d, true);
    print_mc_matrix("mean_FDP [EN+AUG]", snr_labels,
                    rho_labels, en_aug_2d, false);
    print_mc_matrix("mean_TPP [IEN]",    snr_labels,
                    rho_labels, ien_2d,    true);
    print_mc_matrix("mean_FDP [IEN]",    snr_labels,
                    rho_labels, ien_2d,    false);
    save_mc_2d_tables(tag, snr_labels, rho_labels,
                      en_2d, en_aug_2d, ien_2d, cfg, "gvs_" + tag + "_2d");
}


// ==============================================================================
//  main
// ==============================================================================

int main()
{
    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    // ==========================================================================
    //  Simulation parameters
    // ==========================================================================

    GVSSimConfig cfg;
    cfg.n         = 200;
    cfg.p         = 500;
    cfg.sd_x      = 0.1;
    cfg.tFDR      = 0.1;
    cfg.K         = 20;
    cfg.num_MC    = 200;
    cfg.base_seed = 2026;
    cfg.corr_max  = 0.98;
    cfg.snr       = 2.0;    // fixed SNR for Part 2 (rho sweep)

    const std::string scenario_tag = "Neg-Traps";

    // ==========================================================================
    //  Run simulations
    // ==========================================================================

    sim_neg_traps_part1(cfg, scenario_tag);
    sim_neg_traps_part2(cfg, scenario_tag);
    sim_neg_traps_part3(cfg, scenario_tag);

    std::cout << scenario_tag << " GVS MC simulations complete.\n";
    return 0;
}
