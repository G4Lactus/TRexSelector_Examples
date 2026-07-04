// ==============================================================================
// demo_trex_gvs_08_mc_sim_block_bench.cpp
// ==============================================================================
/**
 * @file demo_trex_gvs_08_mc_sim_block_bench.cpp
 *
 * @brief Monte Carlo benchmark: T-Rex+GVS on the block-equicorrelated DGP.
 *
 * @details
 * Tests four method variants over a (rho x snr) grid for each of three
 * truth scenarios:
 *
 *   M1 = T-Rex + EN, HAC-clustered groups
 *   M2 = T-Rex + EN, oracle groups
 *   M3 = T-Rex + IEN, HAC-clustered groups
 *   M4 = T-Rex + IEN, oracle groups
 *
 * Truth scenarios:
 *   - INDIVIDUAL
 *   - REPRESENTATIVE
 *   - WHOLE
 *
 * Fixed benchmark constants:
 *   n = 200, p = 500, G = 100, block_size = 5, n_active_blocks = 10
 *   tFDR = 0.1, K = 20, corr_max = 0.5, b = 3.0
 *   snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}
 *   rho_grid = {0.5, 0.9}
 *   num_MC = 500
 *
 * Seeds are staggered by grid-cell index so that each (rho, snr) cell uses
 * a distinct seed band:
 *
 *   cell_base_seed = base_seed + cell_index * num_MC
 *
 * where cell_index = i_rho * n_snr + i_snr.
 */
// ==============================================================================

// std includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// project includes
#include "trex_gvs_simulation_utils.hpp"

// ==============================================================================

using namespace gvs_demo;

// =============================================================================
// Preset
// =============================================================================

struct BlockBenchPreset {
    int n = 200;
    int p = 500;
    int G = 100;
    int block_size = 5;
    int n_active_blocks = 10;

    double tFDR = 0.1;
    std::size_t K = 20;
    double corr_max = 0.5;
    int base_seed = 2026;
    std::size_t num_MC = 500;
    double b = 3.0;

    std::vector<double> snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
    std::vector<double> rho_grid = {0.5, 0.9};
};

// =============================================================================
// Helpers
// =============================================================================

static BlockBenchPreset make_block_bench_preset()
{
    return BlockBenchPreset{};
}

static std::string scenario_to_tag(BlockScenario scenario)
{
    switch (scenario) {
    case BlockScenario::INDIVIDUAL:
        return "INDIVIDUAL";
    case BlockScenario::REPRESENTATIVE:
        return "REPRESENTATIVE";
    case BlockScenario::WHOLE:
        return "WHOLE";
    }
    return "UNKNOWN";
}

static BlockBenchConfig make_block_bench_cfg(
    const BlockBenchPreset& preset,
    std::size_t i_rho,
    std::size_t i_snr)
{
    BlockBenchConfig cfg;
    cfg.n = preset.n;
    cfg.p = preset.p;
    cfg.G = preset.G;
    cfg.block_size = preset.block_size;
    cfg.n_active_blocks = preset.n_active_blocks;
    cfg.tFDR = preset.tFDR;
    cfg.K = preset.K;
    cfg.num_MC = preset.num_MC;
    cfg.base_seed = preset.base_seed + static_cast<int>(
        (i_rho * preset.snr_grid.size() + i_snr) * preset.num_MC);
    cfg.corr_max = preset.corr_max;
    cfg.b = preset.b;
    return cfg;
}

static BlockBenchConfig make_print_cfg(const BlockBenchPreset& preset)
{
    BlockBenchConfig cfg;
    cfg.n = preset.n;
    cfg.p = preset.p;
    cfg.G = preset.G;
    cfg.block_size = preset.block_size;
    cfg.n_active_blocks = preset.n_active_blocks;
    cfg.tFDR = preset.tFDR;
    cfg.K = preset.K;
    cfg.num_MC = preset.num_MC;
    cfg.base_seed = preset.base_seed;
    cfg.corr_max = preset.corr_max;
    cfg.b = preset.b;
    return cfg;
}

static std::string make_progress_label(
    const std::string& tag,
    double rho,
    double snr,
    std::size_t i_rho,
    std::size_t i_snr,
    std::size_t n_rho,
    std::size_t n_snr)
{
    std::ostringstream oss;
    oss << tag
        << " rho=" << std::fixed << std::setprecision(2) << rho
        << " snr=" << snr
        << " [" << (i_rho * n_snr + i_snr + 1)
        << "/" << (n_rho * n_snr) << "]";
    return oss.str();
}

// =============================================================================
// One cell
// =============================================================================

/**
 * @brief Run one (rho, snr) cell for one truth scenario.
 */
static BlockGridResult run_block_bench_cell(
    const BlockBenchPreset& preset,
    BlockScenario scenario,
    double rho,
    double snr,
    std::size_t i_rho,
    std::size_t i_snr)
{
    const BlockBenchConfig cfg = make_block_bench_cfg(preset, i_rho, i_snr);
    const TRexControlParameter trex_ctrl = make_gvs_trex_control(cfg.K);

    auto dgp_fn = [&](unsigned seed) -> GVSDGPResult {
        return make_block_equicorr_dgp(
            cfg.n,
            cfg.p,
            cfg.G,
            cfg.block_size,
            cfg.n_active_blocks,
            rho,
            snr,
            scenario,
            seed,
            cfg.b);
    };

    const std::string label = make_progress_label(
        scenario_to_tag(scenario),
        rho,
        snr,
        i_rho,
        i_snr,
        preset.rho_grid.size(),
        preset.snr_grid.size());

    return run_block_mc(dgp_fn, cfg, trex_ctrl, label);
}

// =============================================================================
// One scenario
// =============================================================================

/**
 * @brief Run the full (rho x snr) grid for one truth scenario.
 */
static void run_block_bench_scenario(
    const BlockBenchPreset& preset,
    BlockScenario scenario)
{
    const std::string tag = scenario_to_tag(scenario);

    cdiag::print_section_header(
        "Block Bench MC: " + tag +
        "\nn=" + std::to_string(preset.n) +
        " p=" + std::to_string(preset.p) +
        " G=" + std::to_string(preset.G) +
        " block=" + std::to_string(preset.block_size) +
        " active=" + std::to_string(preset.n_active_blocks) +
        " tFDR=" + fmt_num(preset.tFDR) +
        " K=" + std::to_string(preset.K) +
        " MC=" + std::to_string(preset.num_MC) +
        " corr_max=" + fmt_num(preset.corr_max));

    const std::size_t n_rho = preset.rho_grid.size();
    const std::size_t n_snr = preset.snr_grid.size();

    std::vector<std::vector<BlockGridResult>> results(
        n_rho, std::vector<BlockGridResult>(n_snr));

    for (std::size_t i_rho = 0; i_rho < n_rho; ++i_rho) {
        const double rho = preset.rho_grid[i_rho];

        for (std::size_t i_snr = 0; i_snr < n_snr; ++i_snr) {
            const double snr = preset.snr_grid[i_snr];

            results[i_rho][i_snr] =
                run_block_bench_cell(preset, scenario, rho, snr, i_rho, i_snr);
        }
    }

    print_block_grid_table(
        tag,
        scenario,
        preset.rho_grid,
        preset.snr_grid,
        results,
        make_print_cfg(preset));
}

// =============================================================================
// Top-level benchmark
// =============================================================================

/**
 * @brief Run the full block-equicorrelated benchmark for all truth scenarios.
 */
static void run_block_benchmark(const BlockBenchPreset& preset)
{
    cdiag::print_section_header(
        "Block-equicorrelated benchmark"
        "\nDGP = make_block_equicorr_dgp"
        "\nmethods = {EN clustered, EN oracle, IEN clustered, IEN oracle}"
        "\nscenarios = {INDIVIDUAL, REPRESENTATIVE, WHOLE}");

    run_block_bench_scenario(preset, BlockScenario::INDIVIDUAL);
    run_block_bench_scenario(preset, BlockScenario::REPRESENTATIVE);
    run_block_bench_scenario(preset, BlockScenario::WHOLE);

    std::cout << "Block benchmark MC simulations complete.\n";
}

// =============================================================================
// main
// =============================================================================

int main()
{

    // ==========================================================================
    //  Simulation parameters
    // ==========================================================================

    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";


    // ==========================================================================
    //  Run simulations
    // ==========================================================================

    if (true) run_block_benchmark(make_block_bench_preset());

    return 0;
}
