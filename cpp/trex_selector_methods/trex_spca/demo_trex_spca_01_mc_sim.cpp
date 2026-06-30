// ==============================================================================
// demo_trex_spca_01_mc_sim.cpp
// ==============================================================================
/**
 * @file demo_trex_spca_01_mc_sim.cpp
 *
 * @brief Monte Carlo simulation comparing T-Rex SPCA against PCA baselines.
 *
 * @details
 *  DGP: sparse M-factor model — X = Z * V^T + E.
 *       n=50, p=100, M=3 latent factors, p1=5 active loadings per factor,
 *       overlap_pool_size=30 (shared candidate-index pool across factors).
 *
 *  Section 1: SNR sweep   snr_db in {-10, -5, 0, 5, 10} dB, num_MC=100.
 *    Methods compared:
 *      1. OrdPCA             — ordinary PCA (baseline; no sparsity constraint)
 *      2. OraclePCA          — thresholded PCA with known support cardinality p1
 *      3. TRexSPCA-EN-Act    — T-Rex+GVS(EN/TENET, Gram-based) + ActiveSet
 *      4. TRexSPCA-ENaug-Act — T-Rex+GVS(EN/TENETAug, augmented) + ActiveSet
 *      5. TRexSPCA-EN-Thr    — T-Rex+GVS(EN/TENET, Gram-based) + Thresholded
 *      6. TRexSPCA-ENaug-Thr — T-Rex+GVS(EN/TENETAug, augmented) + Thresholded
 *
 *    All methods see the same z-score-standardized X (correlation-PCA footing);
 *    the per-PC T-Rex selector runs with ScalingMode::ZSCORE.
 *
 *    EN augmentation follows lm_dummy.R exactly:
 *      X_aug = d₂·[X; d₁·Iₚ; 0_L],  D_aug = d₂·[D; 0_p; d₁·I_L],
 *      y_aug = [y; 0_{p+L}],  then column-wise centre + L2-normalise.
 *
 *  Metrics (num_MC-averaged; TPR/FDR on PC1 support only — see evaluate_spca):
 *    - FDR  — false discovery rate of PC1's estimated loading support
 *    - TPR  — true positive rate  of PC1's estimated loading support
 *    - PEV  — cumulative percentage of explained variance (all M components)
 */
// ==============================================================================

// std includes
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// OpenMP compatibility layer
#include <utils/openmp/utils_openmp.hpp>

// Demo utilities (transitively includes trex_spca.hpp, pca.hpp, utils_openmp.hpp)
#include "trex_spca_sim_utils.hpp"

// ==============================================================================

// Namespace usage
using namespace spca_sim;

// ==============================================================================
// Monte Carlo Simulation — SNR sweep
// ==============================================================================

/**
 * @brief Run the SNR-sweep MC simulation and save results.
 *
 * @param cfg        Simulation parameters (n, p, p1, M, tFDR, num_MC, seed).
 * @param snr_values SNR grid in dB.
 * @param stem       Output file stem (no folder, no extension).
 */
void demo_trex_spca_mc_snr_sweep(const SimConfig&           cfg,
                                 const std::vector<double>& snr_values,
                                 const std::string&         stem)
{
    std::cout << "================================================================\n";
    std::cout << "  T-Rex SPCA MC Simulation — SNR Sweep\n";
    std::cout << "================================================================\n";
    std::cout << "  n=" << cfg.n << "  p=" << cfg.p << "  M=" << cfg.M
              << "  p1=" << cfg.p1 << "  tFDR=" << cfg.tFDR
              << "  num_MC=" << cfg.num_MC << "\n\n";

    const std::vector<std::string> method_names = {
        "OrdPCA",
        "OraclePCA",
        "TRexSPCA-EN-Act",
        "TRexSPCA-ENaug-Act",
        "TRexSPCA-EN-Thr",
        "TRexSPCA-ENaug-Thr"
    };

    std::map<std::string, Eigen::VectorXd> fdr_map, tpr_map, pev_map;
    for (const auto& name : method_names) {
        fdr_map[name] = Eigen::VectorXd(snr_values.size());
        tpr_map[name] = Eigen::VectorXd(snr_values.size());
        pev_map[name] = Eigen::VectorXd(snr_values.size());
    }

    for (std::size_t si = 0; si < snr_values.size(); ++si) {
        const double   snr     = snr_values[si];
        const int      snr_int = static_cast<int>(std::round(snr));
        const unsigned offset  = static_cast<unsigned>(cfg.base_seed + snr_int);

        std::ostringstream snr_ss;
        snr_ss << std::fixed << std::setprecision(1) << snr;
        const std::string snr_lbl = "SNR=" + snr_ss.str() + " dB";

        std::cout << "--------------------------------------------------\n";
        std::cout << snr_lbl << "\n\n";

        // Thread safe parameter capture all DGP parameters by value for OMP thread safety
        const auto make_data = [cfg, snr](unsigned seed) -> FactorModelData {
            return DataGenerator::generate_sparse_factor_model(
                cfg.n, cfg.p, cfg.p1, cfg.M, snr,
                cfg.overlap_pool_size, seed);
        };

        const auto ei = static_cast<Eigen::Index>(si);

        // 1. Ordinary PCA
        auto r1 = run_mc_trials_pca(
            cfg.num_MC, snr_lbl + " [OrdPCA]", make_data, offset);
        fdr_map["OrdPCA"](ei) = r1.avg_fdr;
        tpr_map["OrdPCA"](ei) = r1.avg_tpr;
        pev_map["OrdPCA"](ei) = r1.avg_pev;


        // 2. Oracle Thresholded PCA
        auto r2 = run_mc_trials_oracle_pca(
            cfg.num_MC, snr_lbl +
            " [OraclePCA]", make_data, cfg.p1, offset
        );

        fdr_map["OraclePCA"](ei) = r2.avg_fdr;
        tpr_map["OraclePCA"](ei) = r2.avg_tpr;
        pev_map["OraclePCA"](ei) = r2.avg_pev;


        // 3. T-Rex SPCA — EN (Gram-based TENET) + ActiveSet
        auto r3 = run_mc_trials_trex_spca(
            cfg.num_MC, snr_lbl + " [TRexSPCA-EN-Act]",
            make_data, cfg.tFDR,
            SPCAMode::ActiveSet,
            offset,
            cfg.lambda_2,
            ENSolverType::TENET,
            ScalingMode::ZSCORE);
        fdr_map["TRexSPCA-EN-Act"](ei) = r3.avg_fdr;
        tpr_map["TRexSPCA-EN-Act"](ei) = r3.avg_tpr;
        pev_map["TRexSPCA-EN-Act"](ei) = r3.avg_pev;


        // 4. T-Rex SPCA — ENaug (augmented TENET) + ActiveSet
        auto r4 = run_mc_trials_trex_spca(
            cfg.num_MC, snr_lbl + " [TRexSPCA-ENaug-Act]",
            make_data, cfg.tFDR,
            SPCAMode::ActiveSet,
            offset,
            cfg.lambda_2,
            ENSolverType::TENET_AUG,
            ScalingMode::ZSCORE);
        fdr_map["TRexSPCA-ENaug-Act"](ei) = r4.avg_fdr;
        tpr_map["TRexSPCA-ENaug-Act"](ei) = r4.avg_tpr;
        pev_map["TRexSPCA-ENaug-Act"](ei) = r4.avg_pev;


        // 5. T-Rex SPCA — EN (Gram-based TENET) + Thresholded
        auto r5 = run_mc_trials_trex_spca(
            cfg.num_MC, snr_lbl + " [TRexSPCA-EN-Thr]",
            make_data, cfg.tFDR,
            SPCAMode::Thresholded,
            offset,
            cfg.lambda_2,
            ENSolverType::TENET,
            ScalingMode::ZSCORE);
        fdr_map["TRexSPCA-EN-Thr"](ei) = r5.avg_fdr;
        tpr_map["TRexSPCA-EN-Thr"](ei) = r5.avg_tpr;
        pev_map["TRexSPCA-EN-Thr"](ei) = r5.avg_pev;


        // 6. T-Rex SPCA — ENaug (augmented TENET) + Thresholded
        auto r6 = run_mc_trials_trex_spca(
            cfg.num_MC, snr_lbl + " [TRexSPCA-ENaug-Thr]",
            make_data, cfg.tFDR,
            SPCAMode::Thresholded,
            offset,
            cfg.lambda_2,
            ENSolverType::TENET_AUG,
            ScalingMode::ZSCORE);
        fdr_map["TRexSPCA-ENaug-Thr"](ei) = r6.avg_fdr;
        tpr_map["TRexSPCA-ENaug-Thr"](ei) = r6.avg_tpr;
        pev_map["TRexSPCA-ENaug-Thr"](ei) = r6.avg_pev;
    }

    save_and_print_spca_mc_results(cfg.num_MC, stem,
                                    snr_values, method_names,
                                    fdr_map, tpr_map, pev_map);
}


// ==============================================================================
// Main
// ==============================================================================

int main() {

    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";


    // =========================================================
    // Section 1: SNR sweep
    // =========================================================
    if (true)
    {
        SimConfig cfg;
        cfg.n                 = 50;
        cfg.p                 = 100;
        cfg.p1                = 5;
        cfg.M                 = 3;
        cfg.overlap_pool_size = 30;
        cfg.tFDR              = 0.10;
        cfg.num_MC            = 200;
        cfg.base_seed         = 42;
        cfg.lambda_2          = -1; // -1 triggers k-fold CV for lambda_2 selection

        const std::vector<double> snr_values = {-10.0}; // TEMP fast-validation

        demo_trex_spca_mc_snr_sweep(cfg, snr_values, "demo_trex_spca_01_mc_sim");
    }

    return 0;
}
