// ==============================================================================
// demo_trex_08_mc_sim_scalability.cpp
// ==============================================================================
/**
 * @file demo_trex_08_mc_sim_scalability.cpp
 *
 * @brief Scalability benchmark for the T-Rex Selector: runtime and memory
 *        usage evaluated in n and p scale, across three base solvers
 *        (TLARS, TOMP, TAFS).
 *
 * @details Runtime and memory benchmark on the fully memory-mapped pipeline:
 *  every MC trial generates X/y directly into mmap backing files
 *  (SyntheticDataMapped over memmap::MemoryMappedMatrix, removed per trial by
 *  an RAII MmapFileGuard), and use_memory_mapping = true additionally
 *  memory-maps the dummy matrices D and serializes solver warm-start state.
 *
 *  Scenarios (see make_benchmark_scenarios()):
 *    A — n =   300, p =   1,000  (baseline)
 *    B — n = 1,000, p =  10,000  (moderate scale)
 *    C — n = 5,000, p = 100,000  (large sample x high-dimensional)
 *
 *  2D sweep: every scenario runs over the SNR grid {0.5, 1.0, 2.0} for the
 *  base solvers TLARS, TOMP, and TAFS (rho_afs = 0.3), with the demo-03 DGP
 *  (variable support of cardinality s = 10 redrawn per trial, unit
 *  coefficients, Normal predictors/noise).
 *
 *  Per grid point (scenario x solver x SNR; num_MC sequential trials so
 *  wall-clock and RSS readings are undistorted): timing breakdown (DGP,
 *  select(), per-T-iteration, total), memory (post-select RSS, process peak
 *  RSS, X/y backing-file sizes), and FDR / TPR / Avg L / Avg T.
 *
 *  Results go to simulation_results/data/ as an aligned table (.txt) and a
 *  tidy CSV, re-saved after each completed scenario — a crash on a larger
 *  scenario preserves the finished ones.
 *
 * @note This simulation will take a substantial amount of time, aovid if possible.
 */
// ==============================================================================

// std includes
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// System includes (memory profiling hooks)
#include <sys/resource.h>
#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

// Eigen includes
#include <Eigen/Dense>

// OpenMP compatibility layer
#include <utils/openmp/utils_openmp.hpp>

// T-Rex Selector includes
#include <trex_selector_methods/trex_core/trex.hpp>
#include <utils/datageneration/utils_datagen.hpp>
#include <utils/eval_metrics/utils_eval_cdiagnostics.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>
#include <utils/memmap/memory_mapped_matrix.hpp>

// Demo utilities
#include "trex_sim_utils.hpp"
using namespace trex_sim;


// ==============================================================================
// Namespace aliases
// ==============================================================================

namespace cdianostics = trex::utils::eval::cdiagnostics;
namespace datagen     = trex::utils::datageneration::datagen;
namespace dummygen    = trex::utils::datageneration::dummygen;
namespace rates       = trex::utils::eval::rates;

// T-Rex types
using trex::trex_selector_methods::trex_core::TRexControlParameter;
using trex::trex_selector_methods::trex_core::TRexSelector;
using trex::trex_selector_methods::trex_core::LLoopStrategy;
using trex::trex_selector_methods::utils::solver_dispatch::SolverTypeForTRex;


// ==============================================================================
// Memory profiling hooks
// ==============================================================================

/**
 * @brief Current resident set size (RSS) of this process in MiB.
 * @return RSS in MiB, or 0.0 if unavailable on this platform.
 */
static double current_rss_mib() {
#if defined(__APPLE__)
    mach_task_basic_info   info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
        return static_cast<double>(info.resident_size) / (1024.0 * 1024.0);
    return 0.0;
#elif defined(__linux__)
    // /proc/self/statm: second field is resident pages
    std::ifstream statm("/proc/self/statm");
    long total_pages = 0, resident_pages = 0;
    if (statm >> total_pages >> resident_pages)
        return static_cast<double>(resident_pages) *
               static_cast<double>(sysconf(_SC_PAGESIZE)) / (1024.0 * 1024.0);
    return 0.0;
#else
    return 0.0;
#endif
}


/**
 * @brief Physical footprint of this process in MiB: anonymous heap including
 *        compressed / swapped-out pages.
 *        Linux approximation: VmRSS + VmSwap from /proc/self/status.
 */
static double current_footprint_mib() {
#if defined(__APPLE__)
    task_vm_info_data_t    info;
    mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_VM_INFO,
                  reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
        return static_cast<double>(info.phys_footprint) / (1024.0 * 1024.0);
    return 0.0;
#elif defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string line;
    double kib = 0.0;
    while (std::getline(status, line)) {
        if (line.rfind("VmRSS:", 0) == 0 || line.rfind("VmSwap:", 0) == 0)
            kib += std::stod(line.substr(line.find(':') + 1));
    }
    return kib / 1024.0;
#else
    return 0.0;
#endif
}


/**
 * @brief Peak resident set size of this process in MiB.
 * @note  Monotone over the process lifetime — for grid points executed later
 *        in the run it reports the maximum seen so far, not a per-point peak.
 */
static double peak_rss_mib() {
    struct rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
#if defined(__APPLE__)
    return static_cast<double>(ru.ru_maxrss) / (1024.0 * 1024.0);  // bytes
#else
    return static_cast<double>(ru.ru_maxrss) / 1024.0;             // KiB
#endif
}


// ==============================================================================
// Scenario descriptors and result types
// ==============================================================================

/** @brief One scalability grid scenario (data size + MC budget).
 *  @note  Every scenario runs the fully memory-mapped pipeline: X/y are
 *         generated into mmap backing files per trial, and D mmap + solver
 *         serialization are enabled via use_memory_mapping = true. */
struct ScalabilityScenario {
    std::string  label;    ///< Scenario tag ("A", "B", "C") used in table/CSV
    Eigen::Index n;        ///< Number of observations
    Eigen::Index p;        ///< Number of predictors
    std::size_t  num_MC;   ///< Number of sequential MC trials
};

/** @brief Measurements from a single MC trial. */
struct ScalabilityTrialResult {
    double dgp_s    = 0.0;   ///< Data generation wall-clock [s]
    double select_s = 0.0;   ///< Selector construction + select() wall-clock [s]
    double fdp      = 0.0;   ///< False discovery proportion
    double tpp      = 0.0;   ///< True positive proportion
    double L        = 0.0;   ///< Calibrated dummy multiplier
    double T        = 0.0;   ///< Stopping time T_stop
    double data_mib = 0.0;   ///< X + y mmap backing-file sizes on disk [MiB]
    double rss_mib  = 0.0;   ///< Current RSS right after select() [MiB]
    double fp_mib   = 0.0;   ///< Physical footprint right after select() [MiB]
};

/** @brief Aggregate result for one grid point (scenario x solver). */
struct ScalabilityGridPointResult {
    double avg_fdr      = 0.0;   ///< Mean false discovery rate
    double avg_tpr      = 0.0;   ///< Mean true positive rate
    double avg_L        = 0.0;   ///< Mean calibrated dummy multiplier
    double avg_T        = 0.0;   ///< Mean stopping time T_stop
    double avg_dgp_s    = 0.0;   ///< Mean data-generation time [s]
    double avg_select_s = 0.0;   ///< Mean select() time (L-loop + T-loop) [s]
    double avg_titer_s  = 0.0;   ///< Mean select time per T-loop iteration [s]
    double avg_total_s  = 0.0;   ///< Mean total trial time [s]
    double max_rss_mib  = 0.0;   ///< Max post-select() RSS across trials [MiB]
    double max_fp_mib   = 0.0;   ///< Max post-select() footprint across trials [MiB]
    double peak_rss_mib = 0.0;   ///< Process peak RSS at grid-point end [MiB]
    double data_mib     = 0.0;   ///< X + y mmap backing-file sizes [MiB]
};


// ==============================================================================
// Scenario lists
// ==============================================================================

/** @brief Full benchmark scenarios.
 *  Each scenario runs the fully memory-mapped pipeline: X/y.
 */
static std::vector<ScalabilityScenario> make_benchmark_scenarios() {
    return {
        {"A", 300,  1'000, /*num_MC=*/10},
        {"B", 1'000, 10'000, /*num_MC=*/10},
        {"C", 5'000,  100'000, /*num_MC=*/10}
    };
}

/** @brief Tiny end-to-end smoke-test scenario (pipeline verification only). */
static std::vector<ScalabilityScenario> make_smoke_scenarios() {
    return {
        {"S", 150, 500, /*num_MC=*/2},
    };
}


// ==============================================================================
// Configuration helpers
// ==============================================================================

/** @brief T-Rex control parameters (base config of demo 02/03) for one grid point. */
static TRexControlParameter make_scalability_trex_ctrl(const DemoSolverInfo& solver)
{
    TRexControlParameter ctrl;
    ctrl.K                        = 20;
    ctrl.max_dummy_multiplier     = 10;
    ctrl.use_max_T_stop           = true;
    ctrl.dummy_distribution       = dummygen::Distribution::Normal();
    ctrl.lloop_strategy           = LLoopStrategy::STANDARD;
    ctrl.tloop_stagnation_stop    = true;
    ctrl.tloop_max_stagnant_steps = 5;
    // Always on: memory-mapped dummy matrices D + solver warm-start state
    // serialized to disk between T-loop iterations. Together with the
    // mmap-generated X/y this makes the whole pipeline memory-mapped.
    ctrl.use_memory_mapping       = true;
    ctrl.solver_type              = solver.solver_type;
    ctrl.solver_params.lambda2    = solver.lambda2;
    ctrl.solver_params.rho_afs    = solver.rho_afs;
    return ctrl;
}

/**
 * @brief Draw a variable active support of given cardinality
 *        (demo 03 pattern: uniform without replacement, redrawn per trial).
 */
static std::vector<std::size_t> draw_support(unsigned seed, Eigen::Index p,
                                             std::size_t cardinality) {
    std::mt19937 rng_sup(seed + 500000u);
    std::vector<std::size_t> pool(static_cast<std::size_t>(p));
    std::iota(pool.begin(), pool.end(), std::size_t{0});
    std::shuffle(pool.begin(), pool.end(), rng_sup);
    pool.resize(cardinality);
    std::sort(pool.begin(), pool.end());
    return pool;
}

/** @brief Unit or N(0,1) coefficients for the active support (demo 03 pattern). */
static std::vector<double> make_coefs(unsigned seed, std::size_t cardinality,
                                      bool rnd_coef) {
    std::vector<double> coefs;
    coefs.reserve(cardinality);
    if (rnd_coef) {
        std::mt19937 rng_coef(seed + 600000u);
        std::normal_distribution<double> nd(0.0, 1.0);
        for (std::size_t i = 0; i < cardinality; ++i) coefs.push_back(nd(rng_coef));
    } else {
        for (std::size_t i = 0; i < cardinality; ++i) coefs.push_back(1.0);
    }
    return coefs;
}


// ==============================================================================
// Single-trial runners (sequential; timing + memory instrumented)
// ==============================================================================

using bench_clock = std::chrono::steady_clock;

static double seconds_since(const bench_clock::time_point& t0) {
    return std::chrono::duration<double>(bench_clock::now() - t0).count();
}

/** @brief One MC trial on the fully memory-mapped pipeline: X/y are generated
 *         directly into mmap backing files (never materialized in RAM) and
 *         fed to the selector through Eigen::Map views.
 *  @note  MmapFileGuard is declared BEFORE SyntheticDataMapped so reverse-LIFO
 *         destruction closes the mmap handles before removing the backing
 *         files — exception-safe (mirrors demo_trex_06/07). */
static ScalabilityTrialResult run_scalability_trial(
    const ScalabilityScenario&      scenario,
    const TRexControlParameter&     trex_ctrl,
    double                          tFDR,
    double                          snr,
    unsigned                        seed,
    const std::vector<std::size_t>& support,
    const std::vector<double>&      coefs,
    const std::string&              mmap_filestem)
{
    ScalabilityTrialResult m;

    const std::string X_path = mmap_filestem + "_X.dat";
    const std::string y_path = mmap_filestem + "_y.dat";

    MmapFileGuard guard{{X_path, y_path}};

    const auto t0 = bench_clock::now();
    datagen::SyntheticDataMapped data(
        X_path, y_path,
        scenario.n, scenario.p,
        support, coefs, snr,
        static_cast<int>(seed));
    m.dgp_s = seconds_since(t0);

    std::error_code ec;
    const auto X_bytes = std::filesystem::file_size(X_path, ec);
    const auto y_bytes = std::filesystem::file_size(y_path, ec);
    if (!ec)
        m.data_mib = static_cast<double>(X_bytes + y_bytes) / (1024.0 * 1024.0);

    // Map views over the mmap backing files — X never fully resides in RAM
    Eigen::Map<Eigen::MatrixXd> X_map(data.getX().data(), data.rows(), data.cols());
    Eigen::Map<Eigen::VectorXd> y_map(data.getY().data(), data.rows());

    const auto t1 = bench_clock::now();
    TRexSelector selector(X_map, y_map, tFDR, trex_ctrl, -1, false);
    selector.select();
    m.select_s = seconds_since(t1);

    m.rss_mib = current_rss_mib();
    m.fp_mib  = current_footprint_mib();

    const auto sel = selector.getSelectedIndices();
    m.fdp = rates::compute_fdp(sel, support);
    m.tpp = rates::compute_tpp(sel, support);
    m.L   = static_cast<double>(selector.getDummyMultiplierL());
    m.T   = static_cast<double>(selector.getTStop());
    return m;
}


// ==============================================================================
// Grid-point loop (scenario x solver): sequential MC trials
// ==============================================================================

static ScalabilityGridPointResult run_scalability_grid_point(
    const ScalabilityScenario& scenario,
    const DemoSolverInfo&      solver,
    double                     tFDR,
    double                     snr,
    unsigned                   base_seed,
    bool                       rnd_coef,
    std::size_t                cardinality_true_support)
{
    const TRexControlParameter trex_ctrl = make_scalability_trex_ctrl(solver);

    ScalabilityGridPointResult agg;
    const std::string mmap_filestem = "demo08_scalability_" + scenario.label;

    std::cout << "  [" << scenario.label << "/" << solver.solver_name
              << "/SNR=" << std::fixed << std::setprecision(1) << snr << "] "
              << "n=" << scenario.n << ", p=" << scenario.p
              << " — running " << scenario.num_MC
              << " sequential MC trials (mmap) ...\n" << std::flush;

    for (std::size_t mc = 0; mc < scenario.num_MC; ++mc) {

        const unsigned trial_seed = base_seed + static_cast<unsigned>(mc);
        const auto support = draw_support(trial_seed, scenario.p,
                                          cardinality_true_support);
        const auto coefs   = make_coefs(trial_seed, cardinality_true_support,
                                        rnd_coef);

        const ScalabilityTrialResult m = run_scalability_trial(
            scenario, trex_ctrl, tFDR, snr, trial_seed,
            support, coefs, mmap_filestem);

        agg.avg_fdr      += m.fdp;
        agg.avg_tpr      += m.tpp;
        agg.avg_L        += m.L;
        agg.avg_T        += m.T;
        agg.avg_dgp_s    += m.dgp_s;
        agg.avg_select_s += m.select_s;
        agg.avg_titer_s  += m.select_s / std::max(m.T, 1.0);
        agg.avg_total_s  += m.dgp_s + m.select_s;
        agg.max_rss_mib   = std::max(agg.max_rss_mib, m.rss_mib);
        agg.max_fp_mib    = std::max(agg.max_fp_mib, m.fp_mib);
        agg.data_mib      = m.data_mib;

        std::cout << "    trial " << (mc + 1) << "/" << scenario.num_MC
                  << ": dgp=" << std::fixed << std::setprecision(2) << m.dgp_s
                  << "s  select=" << m.select_s
                  << "s  FDP=" << std::setprecision(3) << m.fdp
                  << "  TPP=" << m.tpp
                  << "  L=" << std::setprecision(0) << m.L
                  << "  T=" << m.T
                  << "  RSS=" << std::setprecision(1) << m.rss_mib << " MiB"
                  << "  FP=" << m.fp_mib << " MiB\n"
                  << std::flush;
    }

    const double dMC  = static_cast<double>(scenario.num_MC);
    agg.avg_fdr      /= dMC;
    agg.avg_tpr      /= dMC;
    agg.avg_L        /= dMC;
    agg.avg_T        /= dMC;
    agg.avg_dgp_s    /= dMC;
    agg.avg_select_s /= dMC;
    agg.avg_titer_s  /= dMC;
    agg.avg_total_s  /= dMC;
    agg.peak_rss_mib  = peak_rss_mib();

    std::cout << "  [" << scenario.label << "/" << solver.solver_name
              << "/SNR=" << std::fixed << std::setprecision(1) << snr << "] "
              << "done. avg select=" << std::setprecision(2)
              << agg.avg_select_s << "s  FDR=" << std::setprecision(3)
              << agg.avg_fdr << "  TPR=" << agg.avg_tpr << "\n\n" << std::flush;

    return agg;
}


// ==============================================================================
// Output: aligned table (console + .txt) and tidy CSV
// ==============================================================================

/** @brief Metric row descriptor for the results table / CSV. */
struct ScalabilityMetricSpec {
    std::string name;                                  ///< Row / CSV label
    int         precision;                             ///< Table decimals
    double (*get)(const ScalabilityGridPointResult&);  ///< Value extractor
};

/** @brief Metric rows reported per solver (table row order and CSV order). */
static const std::vector<ScalabilityMetricSpec>& scalability_metrics() {
    static const std::vector<ScalabilityMetricSpec> specs = {
        {"FDR",      4, [](const ScalabilityGridPointResult& r) { return r.avg_fdr;      }},
        {"TPR",      4, [](const ScalabilityGridPointResult& r) { return r.avg_tpr;      }},
        {"Avg L",    2, [](const ScalabilityGridPointResult& r) { return r.avg_L;        }},
        {"Avg T",    2, [](const ScalabilityGridPointResult& r) { return r.avg_T;        }},
        {"DGP s",    3, [](const ScalabilityGridPointResult& r) { return r.avg_dgp_s;    }},
        {"Select s", 3, [](const ScalabilityGridPointResult& r) { return r.avg_select_s; }},
        {"s/T-iter", 4, [](const ScalabilityGridPointResult& r) { return r.avg_titer_s;  }},
        {"Total s",  3, [](const ScalabilityGridPointResult& r) { return r.avg_total_s;  }},
        {"RSS MiB",  1, [](const ScalabilityGridPointResult& r) { return r.max_rss_mib;  }},
        {"FootprMiB", 1, [](const ScalabilityGridPointResult& r) { return r.max_fp_mib;  }},
        {"PeakRSS",  1, [](const ScalabilityGridPointResult& r) { return r.peak_rss_mib; }},
        {"X+y MiB",  1, [](const ScalabilityGridPointResult& r) { return r.data_mib;     }},
    };
    return specs;
}

/**
 * @brief Print averaged scalability results as an aligned table (console +
 *        .txt) and write a tidy long-format CSV (demo 02--07 table style;
 *        grid columns are scenarios instead of SNR values).
 */
static void save_and_print_scalability_results(
    const std::string&                      file_stem,
    double                                  tFDR,
    const std::vector<double>&              snr_values,
    const std::vector<ScalabilityScenario>& scenarios,
    const std::vector<DemoSolverInfo>&      solvers,
    const std::map<std::string, std::vector<ScalabilityGridPointResult>>& results)
{
    const std::size_t S = snr_values.size();
    const std::string folder = DEMO_OUTPUT_DIR;
    std::filesystem::create_directories(folder);

    // 1. Open text file (dual output)
    std::ofstream out_file(folder + file_stem + ".txt");
    auto print_dual = [&](const std::string& text) {
        std::cout << text;
        if (out_file.is_open()) out_file << text;
    };

    // 2. Header + scenario info block
    {
        std::ostringstream ss;
        ss << "\n"
           << "======================================================================\n"
           << "=== T-Rex Scalability Results (tFDR = " << std::fixed
           << std::setprecision(2) << tFDR << ", SNR grid = {";
        for (std::size_t i = 0; i < S; ++i)
            ss << (i ? ", " : "") << std::setprecision(1) << snr_values[i];
        ss << "}) ===\n"
           << "======================================================================\n\n";
        ss << "All scenarios: fully memory-mapped pipeline "
              "(X/y mmap + D mmap + solver serialization)\n\n";
        ss << std::left
           << std::setw(10) << "Scenario" << std::setw(9) << "n"
           << std::setw(10) << "p"
           << std::setw(8)  << "num_MC" << "\n"
           << std::string(37, '-') << "\n";
        for (const auto& sc : scenarios) {
            ss << std::left
               << std::setw(10) << sc.label
               << std::setw(9)  << sc.n
               << std::setw(10) << sc.p
               << std::setw(8)  << sc.num_MC << "\n";
        }
        ss << "\n";
        print_dual(ss.str());
    }

    // 3. Table dimensions
    //    Solver name on its own line, metric rows indented beneath it, so the
    //    value columns stay aligned regardless of solver-name length.
    const int indent_w = 2;
    const int metric_w = 23;
    const int col_w    = 14;
    const std::size_t sep_w =
        static_cast<std::size_t>(indent_w + metric_w)
        + static_cast<std::size_t>(col_w) * scenarios.size() * S;

    // 4. Column header + dashed separator (one column per scenario x SNR,
    //    labelled "<scenario>/<snr>")
    {
        std::ostringstream hdr;
        hdr << std::left << std::string(static_cast<std::size_t>(indent_w), ' ')
            << std::setw(metric_w) << "Scenario/SNR";
        for (const auto& sc : scenarios) {
            for (double snr : snr_values) {
                std::ostringstream lbl;
                lbl << sc.label << "/" << std::fixed << std::setprecision(1)
                    << snr;
                hdr << std::right << std::setw(col_w) << lbl.str();
            }
        }
        hdr << "\n" << std::string(sep_w, '-') << "\n";
        print_dual(hdr.str());
    }

    // 5. Data rows: per solver, one row per metric
    for (const auto& solver : solvers) {
        const auto& grid = results.at(solver.solver_name);
        print_dual("\n" + solver.solver_name + "\n");   // name on its own line
        for (const auto& spec : scalability_metrics()) {
            std::ostringstream row;
            row << std::left
                << std::string(static_cast<std::size_t>(indent_w), ' ')
                << std::setw(metric_w) << spec.name;
            // Only the passed scenarios — on incremental saves `grid` also
            // holds default-initialized entries for scenarios not yet run
            for (std::size_t i = 0; i < scenarios.size() * S; ++i)
                row << std::right << std::fixed
                    << std::setprecision(spec.precision)
                    << std::setw(col_w) << spec.get(grid[i]);
            row << "\n";
            print_dual(row.str());
        }
    }
    print_dual("\n");

    // 6. Tidy long-format CSV (scenario, solver, n, p, num_mc, snr, metric, value)
    std::ofstream csv(folder + file_stem + ".csv");
    if (csv.is_open()) {
        csv << "scenario,solver,n,p,num_mc,snr,metric,value\n"
            << std::fixed << std::setprecision(6);
        for (const auto& solver : solvers) {
            const auto& grid = results.at(solver.solver_name);
            for (std::size_t i = 0; i < scenarios.size(); ++i) {
                const auto& sc = scenarios[i];
                for (std::size_t si = 0; si < S; ++si) {
                    for (const auto& spec : scalability_metrics()) {
                        csv << sc.label << "," << solver.solver_name << ","
                            << sc.n << "," << sc.p << ","
                            << sc.num_MC << "," << snr_values[si] << ","
                            << spec.name << ","
                            << spec.get(grid[i * S + si]) << "\n";
                    }
                }
            }
        }
        std::cout << "[Info] CSV results saved to:             "
                  << folder + file_stem + ".csv\n";
    } else {
        std::cout << "[Warning] Could not open CSV file: "
                  << folder + file_stem + ".csv\n";
    }

    // 7. Footer
    if (out_file.is_open()) {
        std::cout << "[Info] Results successfully saved to: "
                  << folder + file_stem + ".txt\n\n";
        out_file.close();
    } else {
        std::cout << "[Warning] Could not open output file: "
                  << folder + file_stem + ".txt\n\n";
    }
}


// ==============================================================================
// Scalability benchmark driver
// ==============================================================================

void demo_TRexSelector_scalability(const std::vector<ScalabilityScenario>& scenarios,
                                   const std::string& file_stem,
                                   bool rnd_coef) {

    std::cout.setf(std::ios::unitbuf);

    cdianostics::print_section_header(
        "Demo: T-Rex Selector Scalability Benchmark  |  TLARS / TOMP / TAFS");

    const std::size_t cardinality_true_support = 10;
    const double tFDR = 0.1;
    const std::vector<double> snr_values = {0.5, 1.0, 2.0};

    // Base solvers to compare (demo 05 configuration)
    const std::vector<DemoSolverInfo> solvers = {
        {SolverTypeForTRex::TLARS, "TLARS"},
        {SolverTypeForTRex::TOMP,  "TOMP"},
        {SolverTypeForTRex::TAFS,  "TAFS", 0.0, 0.3}
    };

    std::cout << "Fully memory-mapped pipeline in all scenarios "
                 "(X/y mmap + D mmap + solver serialization)\n";
    std::cout << "Scenarios: ";
    for (const auto& sc : scenarios)
        std::cout << sc.label << " (n=" << sc.n << ", p=" << sc.p
                  << ", num_MC=" << sc.num_MC << ")  ";
    std::cout << "\nSNR grid = {";
    for (std::size_t i = 0; i < snr_values.size(); ++i)
        std::cout << (i ? ", " : "") << std::fixed << std::setprecision(1)
                  << snr_values[i];
    std::cout << "},  s = " << cardinality_true_support
              << ",  tFDR = " << std::setprecision(2) << tFDR << "\n\n";

    // Results: solver -> grid points flattened as scenario x SNR
    // (index = sc_idx * snr_values.size() + snr_idx)
    const std::size_t S = snr_values.size();
    std::map<std::string, std::vector<ScalabilityGridPointResult>> results;
    for (const auto& solver : solvers)
        results[solver.solver_name].resize(scenarios.size() * S);

    // Scenario-outer loop: every solver finishes the current (smaller)
    // scenario before the next (larger) one starts, and results are re-saved
    // after each completed scenario — an OOM kill or crash on a later
    // scenario preserves everything finished so far.
    for (std::size_t sc_idx = 0; sc_idx < scenarios.size(); ++sc_idx) {
        const auto& scenario = scenarios[sc_idx];

        cdianostics::print_section_header(
            "Scenario " + scenario.label + ": n=" + std::to_string(scenario.n) +
            ", p=" + std::to_string(scenario.p));

        const double X_gib = static_cast<double>(scenario.n) *
                             static_cast<double>(scenario.p) * 8.0 /
                             (1024.0 * 1024.0 * 1024.0);
        std::cout << "  Expected X backing file: " << std::fixed
                  << std::setprecision(2) << X_gib << " GiB "
                  << "(mmap D grows up to 10x that)\n\n" << std::flush;

        // 2D sweep: solver x SNR:
        // The base seed is unique per (scenario, SNR) and shared across solvers
        //  -> every solver sees identical data at each grid point.
        for (const auto& solver : solvers) {
            for (std::size_t snr_idx = 0; snr_idx < S; ++snr_idx) {
                const unsigned base_seed =
                    24u + static_cast<unsigned>(sc_idx * S + snr_idx) * 1000u;
                results[solver.solver_name][sc_idx * S + snr_idx] =
                    run_scalability_grid_point(
                        scenario, solver, tFDR, snr_values[snr_idx],
                        base_seed, rnd_coef, cardinality_true_support);
            }
        }

        // Incremental save over completed scenarios
        const std::vector<ScalabilityScenario> completed(
            scenarios.begin(),
            scenarios.begin() + static_cast<std::ptrdiff_t>(sc_idx) + 1);
        save_and_print_scalability_results(file_stem, tFDR, snr_values,
                                           completed, solvers, results);
    }
    std::cout << "\n\n";
}


// ==============================================================================
// Main
// ==============================================================================

int main() {

    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    // ============================================================
    // Smoke test: tiny sizes, verifies the full pipeline end-to-end
    // ============================================================
    if (false)
        demo_TRexSelector_scalability(
            make_smoke_scenarios(),
            "demo_trex_08_scalability_smoke",
            /*rnd_coef=*/false);

    // ============================================================
    // Scalability benchmark: scenarios A / B / C
    // (sizing rationale in make_benchmark_scenarios())
    // ============================================================
    if (true)
        demo_TRexSelector_scalability(
            make_benchmark_scenarios(),
            "demo_trex_08_scalability_results",
            /*rnd_coef=*/false);

    return 0;
}
