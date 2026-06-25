// ==============================================================================
// demo_trex_04_lloop_strategies.cpp
// ==============================================================================
/**
 * @file demo_trex_04_lloop_strategies.cpp
 *
 * @brief Monte Carlo simulation comparing the six L-loop strategies of the
 *        T-Rex Selector, with TLARS held fixed as the base solver.
 *
 * @details Mirrors R/trex_selector_methods/trex/demo_trex_04_lloop_strategies.R
 *
 *  Runs in two modes controlled by the `block_support` parameter:
 *    block_support = true  — contiguous block support {0, 1, ..., s-1}.
 *    block_support = false — scattered support redrawn from uniform per trial.
 *
 *  Holds the solver fixed (TLARS) and sweeps over the six L-loop strategies:
 *
 *  Demo content:
 *    STANDARD           — Fresh i.i.d. dummy matrix at each L-loop iteration.
 *    HCONCAT            — Horizontally expand dummy columns.
 *    PERMUTATION        — Re-use the base dummy matrix via column permutations.
 *    PERMUTATION_DIRECT — Seed-based permutations; no base matrix in memory.
 *    DIRECT             — Seed-based i.i.d. draws; no base matrix in memory.
 *    SKIPL              — Skip the L-loop entirely and uses L = max_dummy_multiplier.
 *
 *    High-dimensional  (n = 300, p = 1000, s = 10).
 *    SNR sweep: {0.1, 0.2, ..., 2.0, 5.0}  (21 values).
 *    Two support scenarios:
 *      random support  — redrawn per trial;
 *      block support   — contiguous block {0, 1, ..., s-1}.
 *    Reports averaged FDR / TPR / Avg L / Avg T per L-strategy x SNR.
 */
// ==============================================================================

// std includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// OpenMP compatibility layer
#include <utils/openmp/utils_openmp.hpp>

// T-Rex Selector includes
#include <trex_selector_methods/trex_core/trex.hpp>
#include <utils/datageneration/utils_datagen.hpp>
#include <utils/eval_metrics/utils_eval_cdiagnostics.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>

// Demo utilities
#include "trex_sim_utils.hpp"
using namespace trex_sim;


// ==============================================================================
// Namespace aliases
// ==============================================================================

namespace cdianostics = trex::utils::eval::cdiagnostics;
namespace datagen = trex::utils::datageneration::datagen;
namespace dummygen = trex::utils::datageneration::dummygen;

// T-Rex types
using trex::trex_selector_methods::trex_core::TRexControlParameter;
using trex::trex_selector_methods::trex_core::LLoopStrategy;
using trex::trex_selector_methods::utils::solver_dispatch::SolverTypeForTRex;


// ==============================================================================
// L-loop strategy descriptor
// ==============================================================================

/**
 * @brief Descriptor for one L-loop strategy variant.
 */
struct LLoopStrategyInfo {
    /** @brief Row label used in output table and CSV. */
    std::string name;
    /** @brief Enum value passed to TRexControlParameter. */
    LLoopStrategy strategy;
    /**
     * @brief Upper bound for the adaptive L-loop (non-SKIPL) or the fixed L
     *        value for SKIPL (expressed as a multiple of p, so total dummies
     *        = max_dummy_multiplier * p).
     */
    std::size_t max_dummy_multiplier = 10;
};


/** @brief Returns the canonical list of L-loop strategy variants to compare.
 *
 *  Adaptive strategies (STANDARD … DIRECT) all use max_dummy_multiplier = 10 (default).
 *  SKIPL is compared at four explicit L levels: 5p, 10p, 20p, and 50p.
 */
static std::vector<LLoopStrategyInfo> make_lloop_strategies() {
    return {
        //{"STANDARD",           LLoopStrategy::STANDARD},
        //{"HCONCAT",            LLoopStrategy::HCONCAT},
        //{"PERMUTATION",        LLoopStrategy::PERMUTATION},
        //{"PERMUTATION_DIRECT", LLoopStrategy::PERMUTATION_DIRECT},
        //{"DIRECT",             LLoopStrategy::DIRECT},
        {"SKIPL_5p",           LLoopStrategy::SKIPL, 5},
        {"SKIPL_10p",          LLoopStrategy::SKIPL, 10},
        {"SKIPL_20p",          LLoopStrategy::SKIPL, 20}
    };
}


// ==============================================================================
// Monte Carlo Simulation — L-loop strategy comparison (TLARS fixed)
// ==============================================================================

void demo_TRexSelector_lloop_strategies(std::size_t num_MC, bool rnd_coef,
                                        bool block_support) {

    std::cout.setf(std::ios::unitbuf);

    cdianostics::print_section_header(
        "Demo: T-Rex Selector MC \u2014 L-Strategy Comparison  |  TLARS");

    // -----------------------------------------------------------------------
    // Scenario parameters (always high-dimensional, mirrors R file)
    // -----------------------------------------------------------------------
    const Eigen::Index n = 300;
    const Eigen::Index p = 1000;
    const std::size_t  s = 10;
    const double       tFDR = 0.1;

    std::cout << "High-dimensional (p > n)\n";
    std::cout << "n = " << n << ",  p = " << p
              << ",  s = " << s << ",  num_MC = " << num_MC
              << "  [" << (block_support ? "block" : "random") << " support]\n\n";

    // SNR values: 0.1, 0.2, ..., 2.0, 5.0  (21 values)
    std::vector<double> snr_values(20);
    std::iota(snr_values.begin(), snr_values.end(), 1);
    for (auto& x : snr_values) x *= 0.1;
    snr_values.push_back(5.0);

    // -----------------------------------------------------------------------
    // L-loop strategies to sweep
    // -----------------------------------------------------------------------
    const std::vector<LLoopStrategyInfo> strategies = make_lloop_strategies();

    // -----------------------------------------------------------------------
    // Result containers: strategy x SNR
    // -----------------------------------------------------------------------
    std::map<std::string, Eigen::VectorXd> fdr_results_map;
    std::map<std::string, Eigen::VectorXd> tpr_results_map;
    std::map<std::string, Eigen::VectorXd> average_L_results_map;
    std::map<std::string, Eigen::VectorXd> average_T_results_map;

    std::vector<std::string> strategy_names;
    strategy_names.reserve(strategies.size());
    for (const auto& strat : strategies) {
        strategy_names.push_back(strat.name);
        fdr_results_map[strat.name]       = Eigen::VectorXd(snr_values.size());
        tpr_results_map[strat.name]       = Eigen::VectorXd(snr_values.size());
        average_L_results_map[strat.name] = Eigen::VectorXd(snr_values.size());
        average_T_results_map[strat.name] = Eigen::VectorXd(snr_values.size());
    }

    // -----------------------------------------------------------------------
    // Strategy loop
    // -----------------------------------------------------------------------
    for (const auto& strat : strategies) {

        cdianostics::print_section_header("L-strategy: " + strat.name);

        TRexControlParameter ctrl;
        ctrl.solver_type          = SolverTypeForTRex::TLARS;
        ctrl.K                    = 20;
        ctrl.max_dummy_multiplier = strat.max_dummy_multiplier;
        ctrl.use_max_T_stop       = true;
        ctrl.dummy_distribution   = dummygen::Distribution::Normal();
        ctrl.lloop_strategy       = strat.strategy;

        // Simulation loop: sweep over SNR values
        for (std::size_t snr_idx = 0; snr_idx < snr_values.size(); ++snr_idx) {

            const double   snr       = snr_values[snr_idx];
            const unsigned base_seed = 24u + static_cast<unsigned>(snr_idx) * 1000u;

            std::cout << "SNR = " << std::fixed << std::setprecision(2) << snr << "\n";

            const auto make_data = [=](unsigned seed) -> trex_sim::TrexDGPData {
                std::vector<std::size_t> pool;
                if (block_support) {
                    pool.resize(s);
                    std::iota(pool.begin(), pool.end(), std::size_t{0});
                } else {
                    std::mt19937 rng_sup(seed + 500000u);
                    pool.resize(static_cast<std::size_t>(p));
                    std::iota(pool.begin(), pool.end(), std::size_t{0});
                    std::shuffle(pool.begin(), pool.end(), rng_sup);
                    pool.resize(s);
                    std::sort(pool.begin(), pool.end());
                }

                std::vector<double> coefs;
                coefs.reserve(s);
                if (rnd_coef) {
                    std::mt19937 rng_coef(seed + 600000u);
                    std::normal_distribution<double> nd(0.0, 1.0);
                    for (std::size_t i = 0; i < s; ++i)
                        coefs.push_back(nd(rng_coef));
                } else {
                    for (std::size_t i = 0; i < s; ++i)
                        coefs.push_back(1.0);
                }

                datagen::SyntheticData data(
                    n, p,
                    pool, coefs, snr,
                    static_cast<int>(seed),
                    datagen::predictor_policy::Normal(),
                    datagen::noisegen::noise_policy::Normal());

                return {data.getX(), data.getY(), pool};
            };

            std::ostringstream lbl;
            lbl << "SNR=" << std::fixed << std::setprecision(2) << snr
                << " [" << strat.name << "]";
            const auto result = trex_sim::run_mc_trials_trex(
                num_MC, lbl.str(), make_data, tFDR, ctrl,
                base_seed);

            const auto ei = static_cast<Eigen::Index>(snr_idx);
            fdr_results_map[strat.name](ei)        = result.avg_fdr;
            tpr_results_map[strat.name](ei)        = result.avg_tpr;
            average_L_results_map[strat.name](ei)  = result.avg_L;
            average_T_results_map[strat.name](ei)  = result.avg_T;
        }

        std::cout << "\n";
    }

    // -----------------------------------------------------------------------
    // Print and save results
    // -----------------------------------------------------------------------
    const std::string stem = "demo_trex_04_lloop_strategies_results_n"
        + std::to_string(static_cast<std::size_t>(n))
        + "_p" + std::to_string(static_cast<std::size_t>(p))
        + (block_support ? "_block_support" : "_random_support");

    save_and_print_mc_results(
        num_MC,
        stem,
        snr_values,
        strategy_names,
        fdr_results_map,
        tpr_results_map,
        average_L_results_map,
        average_T_results_map);

    std::cout << "\n\n";
}


// ==============================================================================
// Main
// ==============================================================================

int main() {

    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    if (true)
        demo_TRexSelector_lloop_strategies(/*num_MC=*/10,
                                           /*rnd_coef=*/false,
                                           /*block_support=*/false);
    if (false)
        demo_TRexSelector_lloop_strategies(/*num_MC=*/10,
                                           /*rnd_coef=*/false,
                                           /*block_support=*/true);
    return 0;
}
