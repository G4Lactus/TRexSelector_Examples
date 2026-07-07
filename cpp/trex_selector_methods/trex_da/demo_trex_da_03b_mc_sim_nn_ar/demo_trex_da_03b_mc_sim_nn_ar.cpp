// ==============================================================================
// demo_trex_da_03b_mc_sim_nn_ar.cpp
// ==============================================================================
/**
 * @file demo_trex_da_03b_mc_sim_nn_ar.cpp
 *
 * @brief DA+NN selector applied to AR(1) data — Monte Carlo simulations.
 *
 * @details
 *  Uses an AR(1) DGP but runs DA-TRex with DAMethod::NN, exploring whether the
 *  nearest-neighbour aggregation can track the AR(1) correlation structure.
 *  Companion to demo_trex_da_03_mc_sim_nn.cpp (NN DGP).
 *
 *  Part 1: MC SNR sweep (n=300, p=1000, rho=0.7; CappedSpread and Random
 *          support; 3 solvers + base T-Rex comparison).
 *  Part 2: MC rho sweep (n=300, p=1000, SNR=2.0; dual support).
 *  Part 3: 2D SNR × rho sweep (n=300, p=1000; all 3 solvers; CappedSpread
 *          and Random sub-sections).
 *
 * Results are saved to simulation_results/ as TXT and CSV files.
 */
// ==============================================================================

#include "trex_da_sim_common.hpp"

// ==============================================================================

using namespace da_sim;

// Shared MC parameters
namespace {
    constexpr int    mc_n      = 300;
    constexpr int    mc_p      = 1000;
    constexpr int    mc_s      = 10;
    constexpr double mc_rho    = 0.7;
    constexpr double mc_snr    = 2.0;
    constexpr double tFDR      = 0.2;
    constexpr int    K         = 20;
    constexpr int    num_MC    = 200;
    constexpr int    base_seed = 2026;
    constexpr int    max_gap   = 20;
} // namespace


// ==============================================================================
//  Part 1: MC SNR sweep — AR(1) data with DA+NN
//  Fixed: rho=0.7.
//  CappedSpread(max_gap=20) and Random support.
// ==============================================================================

void demo_nn_ar_mc_snr_sweep()
{
    const std::vector<double> snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::NN;

    const auto support_fixed = capped_spread_support(mc_s, mc_p, max_gap);

    const std::string hdr_base =
        "AR(1)+NN rho=" + fmt_num(mc_rho)
        + ", n=" + std::to_string(mc_n)
        + ", p=" + std::to_string(mc_p);

    // A — CappedSpread
    run_snr_sweep(
        "da_nn_ar_snr_capped", snr_grid,
        num_MC, tFDR, base_seed, default_solvers(),
        [&](double snr, unsigned seed) {
            return dgp_ar1(mc_n, mc_p, support_fixed,
                           /*amplitude=*/3.0, mc_rho, snr, seed);
        },
        da_ctrl, make_trex_control(K),
        hdr_base + ", support=CappedSpread(max_gap=" +
            std::to_string(max_gap) + ")",
        /*include_base_trex=*/true);

    // B — Random (redrawn each trial)
    run_snr_sweep(
        "da_nn_ar_snr_random", snr_grid,
        num_MC, tFDR, base_seed, default_solvers(),
        [&](double snr, unsigned seed) {
            return dgp_ar1(mc_n, mc_p,
                           make_support(SupportPolicy::Random,
                                        mc_s, mc_p, 0, seed),
                           /*amplitude=*/3.0, mc_rho, snr, seed);
        },
        da_ctrl, make_trex_control(K),
        hdr_base + ", support=Random (per-trial)",
        /*include_base_trex=*/true);
}


// ==============================================================================
//  Part 2: MC rho sweep — AR(1) data with DA+NN
//  Fixed: SNR=2.0.
//  CappedSpread(max_gap=20) and Random support.
// ==============================================================================

void demo_nn_ar_mc_rho_sweep()
{
    const std::vector<double> rho_grid = {
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::NN;

    const auto support_fixed = capped_spread_support(mc_s, mc_p, max_gap);

    const std::string hdr_base =
        "AR(1)+NN SNR=" + fmt_num(mc_snr)
        + ", n=" + std::to_string(mc_n)
        + ", p=" + std::to_string(mc_p);

    // A — CappedSpread
    run_param_sweep(
        "da_nn_ar_rho_capped", "Rho", rho_grid,
        num_MC, tFDR, base_seed, default_solvers(),
        [&](double rho, unsigned seed) {
            return dgp_ar1(mc_n, mc_p, support_fixed,
                           /*amplitude=*/3.0, rho, mc_snr, seed);
        },
        da_ctrl, make_trex_control(K),
        hdr_base + ", support=CappedSpread(max_gap=" +
            std::to_string(max_gap) + ")",
        /*include_base_trex=*/true);

    // B — Random
    run_param_sweep(
        "da_nn_ar_rho_random", "Rho", rho_grid,
        num_MC, tFDR, base_seed, default_solvers(),
        [&](double rho, unsigned seed) {
            return dgp_ar1(mc_n, mc_p,
                           make_support(SupportPolicy::Random,
                                        mc_s, mc_p, 0, seed),
                           /*amplitude=*/3.0, rho, mc_snr, seed);
        },
        da_ctrl, make_trex_control(K),
        hdr_base + ", support=Random (per-trial)",
        /*include_base_trex=*/true);
}


// ==============================================================================
//  Part 3: 2D SNR × rho sweep — AR(1) data with DA+NN
//  Fixed: amplitude=3.0.
//  Outer loop: solver (TLARS, TAFS, TOMP).
//  Section A: CappedSpread(max_gap=20);
//  Section B: Random (per-trial).
// ==============================================================================

void demo_nn_ar_mc_snr_rho_sweep()
{
    const std::vector<double> snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
    const std::vector<double> rho_grid = {
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    const auto n_snr = snr_grid.size();
    const auto n_rho = rho_grid.size();
    const auto n_s_idx = static_cast<Eigen::Index>(n_snr);
    const auto n_r_idx = static_cast<Eigen::Index>(n_rho);

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::NN;

    const auto support_cs = capped_spread_support(mc_s, mc_p, max_gap);
    const auto solvers    = default_solvers();

    // Per-solver result maps [n_snr × n_rho]
    std::map<std::string, Eigen::MatrixXd> tpp_cs_map,    fdp_cs_map,
                                           sd_tpp_cs_map,  sd_fdp_cs_map;
    std::map<std::string, Eigen::MatrixXd> tpp_rand_map,   fdp_rand_map,
                                           sd_tpp_rand_map, sd_fdp_rand_map;

    for (const auto& sv : solvers) {
        tpp_cs_map[sv.name]     = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
        fdp_cs_map[sv.name]     = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
        sd_tpp_cs_map[sv.name]  = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
        sd_fdp_cs_map[sv.name]  = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
        tpp_rand_map[sv.name]    = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
        fdp_rand_map[sv.name]    = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
        sd_tpp_rand_map[sv.name] = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
        sd_fdp_rand_map[sv.name] = Eigen::MatrixXd::Zero(n_s_idx, n_r_idx);
    }

    // --- Header ---
    cdianostics::print_section_header("AR(1)+NN Part 3: 2D SNR × Rho sweep");
    std::cout << "  n=" << mc_n << "  p=" << mc_p << "  s=" << mc_s
              << "  MC=" << num_MC
              << "  tFDR=" << tFDR << "  amplitude=3.0\n";
    std::cout << "  snr_grid: {";
    for (std::size_t i = 0; i < n_snr; ++i)
        std::cout << (i ? ", " : "")
                  << std::fixed << std::setprecision(2) << snr_grid[i];
    std::cout << "}\n  rho_grid: {";
    for (std::size_t i = 0; i < n_rho; ++i)
        std::cout << (i ? ", " : "")
                  << std::fixed << std::setprecision(1) << rho_grid[i];
    std::cout << "}\n\n";

    // ---- Loop over solvers ----
    for (const auto& sv : solvers) {
        auto trex_ctrl = make_trex_control(K);
        trex_ctrl.solver_type           = sv.type;
        trex_ctrl.solver_params.rho_afs = sv.rho_afs;
        trex_ctrl.tloop_stagnation_stop = sv.use_stagnation;

        std::cout << "\n  Solver: " << sv.name << "\n\n";

        const int total_cells = static_cast<int>(n_snr * n_rho);

        // ---- Section A: CappedSpread(max_gap=20) ----
        std::cout << "  [A] CappedSpread(max_gap=" << max_gap << ") ...\n\n";
        int cell_idx = 0;
        for (std::size_t is = 0; is < n_snr; ++is) {
            const double snr_val = snr_grid[is];
            for (std::size_t ir = 0; ir < n_rho; ++ir) {
                ++cell_idx;
                const double rho_val = rho_grid[ir];

                std::ostringstream lbl;
                lbl << "[" << cell_idx << "/" << total_cells << "] "
                    << sv.name
                    << "  SNR=" << std::fixed << std::setprecision(2) << snr_val
                    << "  rho=" << std::fixed << std::setprecision(1) << rho_val;

                auto res = run_mc_trials(
                    num_MC, lbl.str(),
                    [&](unsigned seed) {
                        return dgp_ar1(mc_n, mc_p, support_cs,
                                       /*amplitude=*/3.0,
                                       rho_val, snr_val, seed);
                    },
                    tFDR, trex_ctrl, da_ctrl,
                    static_cast<unsigned>(base_seed));

                const auto si = static_cast<Eigen::Index>(is);
                const auto ri = static_cast<Eigen::Index>(ir);
                tpp_cs_map[sv.name](si, ri)    = res.avg_tpr;
                fdp_cs_map[sv.name](si, ri)    = res.avg_fdr;
                sd_tpp_cs_map[sv.name](si, ri) = res.sd_tpr;
                sd_fdp_cs_map[sv.name](si, ri) = res.sd_fdr;
            }
        }

        // ---- Section B: Random support (redrawn each trial) ----
        std::cout << "\n  [B] Random support ...\n\n";
        cell_idx = 0;
        for (std::size_t is = 0; is < n_snr; ++is) {
            const double snr_val = snr_grid[is];
            for (std::size_t ir = 0; ir < n_rho; ++ir) {
                ++cell_idx;
                const double rho_val = rho_grid[ir];

                std::ostringstream lbl;
                lbl << "[" << cell_idx << "/" << total_cells << "] "
                    << sv.name
                    << "  SNR=" << std::fixed << std::setprecision(2) << snr_val
                    << "  rho=" << std::fixed << std::setprecision(1) << rho_val
                    << "  Random";

                auto res = run_mc_trials(
                    num_MC, lbl.str(),
                    [&](unsigned seed) {
                        return dgp_ar1(mc_n, mc_p,
                                       make_support(SupportPolicy::Random,
                                                    mc_s, mc_p, 0, seed),
                                       /*amplitude=*/3.0,
                                       rho_val, snr_val, seed);
                    },
                    tFDR, trex_ctrl, da_ctrl,
                    static_cast<unsigned>(base_seed));

                const auto si = static_cast<Eigen::Index>(is);
                const auto ri = static_cast<Eigen::Index>(ir);
                tpp_rand_map[sv.name](si, ri)    = res.avg_tpr;
                fdp_rand_map[sv.name](si, ri)    = res.avg_fdr;
                sd_tpp_rand_map[sv.name](si, ri) = res.sd_tpr;
                sd_fdp_rand_map[sv.name](si, ri) = res.sd_fdr;
            }
        }
    }

    // ---- Save results ----
    const std::string header_extra =
        "n=" + std::to_string(mc_n) +
        ", p=" + std::to_string(mc_p) +
        ", s=" + std::to_string(mc_s) +
        ", tFDR=" + fmt_num(tFDR) +
        ", amplitude=3.0";

    save_and_print_2d_two_support_results(
        "da_nn_ar_snr_rho",
        "SNR", "Rho",
        snr_grid, rho_grid,
        num_MC,
        solvers,
        tpp_cs_map,   fdp_cs_map,   sd_tpp_cs_map,  sd_fdp_cs_map,
        tpp_rand_map, fdp_rand_map, sd_tpp_rand_map, sd_fdp_rand_map,
        header_extra);
}


// ==============================================================================
//  main
// ==============================================================================

int main() {

    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    if (true) demo_nn_ar_mc_snr_sweep();
    if (true) demo_nn_ar_mc_rho_sweep();
    if (true) demo_nn_ar_mc_snr_rho_sweep();

    std::cout << "\nAR(1)+NN MC simulation complete.\n";
    return 0;
}
