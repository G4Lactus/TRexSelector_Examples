// ==============================================================================
// demo_trex_da_02_mc_sim_equi_and_bt.cpp
// ==============================================================================
/**
 * @file demo_trex_da_02_mc_sim_equi_and_bt.cpp
 *
 * @brief DA-TRex Monte Carlo simulations for equicorrelated and BT DGPs.
 *
 * @details
 *  Covers the R references demo_trex_da_02_equi.R and
 *  demo_trex_da_03_bt_equi_blocks.R.
 *
 *  Part 1: MC SNR sweep — Equi DGP (n=300, p=1000, rho=0.25, Random
 *          support, 3 solvers + base T-Rex comparison).
 *
 *  Part 2: MC rho sweep — Equi DGP (n=300, p=1000, SNR=2.0, Random
 *          support, 3 solvers + base T-Rex comparison).
 *
 *  Part 3: MC SNR sweep — BT DGP (n=300, p=1000, n_blocks=10,
 *          rho_within=0.5, rho_between=0.1, linkage=Single, OnePerBlock
 *          support, 3 solvers + base T-Rex comparison).
 *
 *  Part 4: MC linkage sweep — BT DGP (n=300, p=1000, n_blocks=10,
 *          rho_within=0.8, rho_between=0.1, OnePerBlock support; outer
 *          loop over Single / Complete / Average linkage).
 *
 * Results are saved to simulation_results/ as TXT and CSV files.
 */
// ==============================================================================

#include "trex_da_sim_common.hpp"

// ==============================================================================

using namespace da_sim;

// ==============================================================================
//  Part 1: MC SNR sweep — Equicorrelated DGP
// ==============================================================================

void demo_equi_mc_snr_sweep()
{
    constexpr int    n         = 300;
    constexpr int    p         = 1000;
    constexpr int    s         = 10;
    constexpr int    K         = 20;
    constexpr int    num_MC    = 200;
    constexpr double tFDR      = 0.2;
    constexpr int    base_seed = 2026;
    constexpr double rho       = 0.25;
    const auto       policy = SupportPolicy::Random;

    const std::vector<double> snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::EQUI;

    run_snr_sweep(
        "da_equi_snr",
        snr_grid,
        num_MC,
        tFDR,
        base_seed,
        default_solvers(),
        [&](double snr, unsigned seed) {
            return dgp_equi(n, p,
                            make_support(policy, s, p, 0, seed),
                            /*amplitude=*/3.0,
                            rho, snr, seed);
        },
        da_ctrl,
        make_trex_control(K),
        "Equi rho=" + fmt_num(rho)
        + ", n=" + std::to_string(n)
        + ", p=" + std::to_string(p)
        + ", support=" + support_policy_label(policy),
        /*include_base_trex=*/true);
}


// ==============================================================================
//  Part 2: MC rho sweep — Equicorrelated DGP
// ==============================================================================

void demo_equi_mc_rho_sweep()
{
    constexpr int    n         = 300;
    constexpr int    p         = 1000;
    constexpr int    s         = 10;
    constexpr int    K         = 20;
    constexpr int    num_MC    = 200;
    constexpr double tFDR      = 0.2;
    constexpr int    base_seed = 2026;
    constexpr double snr       = 2.0;
    const auto       policy   = SupportPolicy::Random;

    const std::vector<double> rho_grid = {
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::EQUI;

    run_param_sweep(
        "da_equi_rho",
        "Rho",
        rho_grid,
        num_MC,
        tFDR,
        base_seed,
        default_solvers(),
        [&](double rho, unsigned seed) {
            return dgp_equi(n, p,
                            make_support(policy, s, p, 0, seed),
                            /*amplitude=*/3.0,
                            rho, snr, seed);
        },
        da_ctrl,
        make_trex_control(K),
        "Equi SNR=" + fmt_num(snr)
        + ", n=" + std::to_string(n)
        + ", p=" + std::to_string(p)
        + ", support=" + support_policy_label(policy),
        /*include_base_trex=*/true);
}


// ==============================================================================
//  Part 3: MC SNR sweep — BT (hierarchical-block) DGP, linkage = Single
// ==============================================================================

void demo_bt_mc_snr_sweep()
{
    constexpr int    n            = 300;
    constexpr int    p            = 1000;
    constexpr int    s            = 10;
    constexpr int    K            = 20;
    constexpr int    num_MC       = 200;
    constexpr double tFDR         = 0.2;
    constexpr int    base_seed    = 2026;
    constexpr int    n_blocks     = 10;
    constexpr double rho_within   = 0.5;
    constexpr double rho_between  = 0.1;
    const auto       policy       = SupportPolicy::OnePerBlock;

    const std::vector<double> snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    // Fixed support (OnePerBlock; seed-independent for fixed s/p/n_blocks)
    const auto support_fixed = make_support(policy, s, p, 0, 0u);

    auto da_ctrl = make_da_bt_control();
    da_ctrl.hc_linkage = hac::LinkageMethod::Single;

    // Stem is deliberately distinct from Part 4's per-linkage stems
    // ("da_bt_snr_Single" etc.): on a case-insensitive filesystem
    // "da_bt_snr_single" and "da_bt_snr_Single" would otherwise be the same
    // file and silently overwrite each other.
    run_snr_sweep(
        "da_bt_snr_single_baseline",
        snr_grid,
        num_MC,
        tFDR,
        base_seed,
        default_solvers(),
        [&](double snr, unsigned seed) {
            return dgp_bt(n, p, support_fixed,
                          /*amplitude=*/3.0,
                          n_blocks, rho_within, rho_between,
                          snr, seed);
        },
        da_ctrl,
        make_trex_control(K),
        "BT n_blocks=" + std::to_string(n_blocks)
        + " rho_w=" + fmt_num(rho_within)
        + " rho_b=" + fmt_num(rho_between)
        + ", n=" + std::to_string(n) + ", p=" + std::to_string(p)
        + ", linkage=Single"
        + ", support=" + support_policy_label(policy),
        /*include_base_trex=*/true);
}


// ==============================================================================
//  Part 4: MC linkage sweep — BT DGP, high within-correlation
//  Outer loop: Single / Complete / Average linkage methods
// ==============================================================================

void demo_bt_mc_linkage_sweep()
{
    constexpr int    n           = 300;
    constexpr int    p           = 1000;
    constexpr int    s           = 10;
    constexpr int    K           = 20;
    constexpr int    num_MC      = 200;
    constexpr double tFDR        = 0.2;
    constexpr int    base_seed   = 2026;
    constexpr int    n_blocks    = 10;
    constexpr double rho_within  = 0.8;
    constexpr double rho_between = 0.1;
    const auto       policy      = SupportPolicy::OnePerBlock;

    const std::vector<double> snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    const auto support_fixed = make_support(policy, s, p, 0, 0u);

    const std::vector<std::pair<std::string, hac::LinkageMethod>> linkages = {
        {"Single",   hac::LinkageMethod::Single},
        {"Complete", hac::LinkageMethod::Complete},
        {"Average",  hac::LinkageMethod::Average}
    };

    for (const auto& [lnk_str, lnk_val] : linkages) {
        auto da_ctrl = make_da_bt_control();
        da_ctrl.hc_linkage = lnk_val;

        run_snr_sweep(
            "da_bt_snr_" + lnk_str,
            snr_grid,
            num_MC,
            tFDR,
            base_seed,
            default_solvers(),
            [&](double snr, unsigned seed) {
                return dgp_bt(n, p, support_fixed,
                              /*amplitude=*/3.0,
                              n_blocks, rho_within, rho_between,
                              snr, seed);
            },
            da_ctrl,
            make_trex_control(K),
            "BT n_blocks=" + std::to_string(n_blocks)
            + " rho_w=" + fmt_num(rho_within)
            + " rho_b=" + fmt_num(rho_between)
            + ", n=" + std::to_string(n) + ", p=" + std::to_string(p)
            + ", linkage=" + lnk_str
            + ", support=" + support_policy_label(policy),
            /*include_base_trex=*/false);
    }
}


// ==============================================================================
//  main
// ==============================================================================

int main() {

    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    if (true) demo_equi_mc_snr_sweep();
    if (true) demo_equi_mc_rho_sweep();
    if (true) demo_bt_mc_snr_sweep();
    if (true) demo_bt_mc_linkage_sweep();

    std::cout << "\nEqui + BT MC simulation complete.\n";
    return 0;
}
