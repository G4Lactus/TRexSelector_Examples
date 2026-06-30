// ==============================================================================
// validation_trex_spca_01_lambda2_probe.cpp
// ==============================================================================
/**
 * @file validation_trex_spca_01_lambda2_probe.cpp
 *
 * @brief Isolated diagnostic for the lambda_2 (ridge penalty) determination
 *        used by T-Rex SPCA, decoupled from the full selection pipeline.
 *
 * @details
 *  The full SPCA Monte Carlo (demo_01) hides lambda_2 inside the selector.
 *  This probe reproduces *exactly* what TRexGVSSelector::computeLambda2() does
 *  for PC1, and prints every intermediate quantity so it can be compared
 *  number-for-number against R's glmnet::cv.glmnet(alpha = 0).
 *
 *  Production path being reproduced (see trex_spca.cpp / trex_gvs.cpp):
 *    1. X is centered column-wise (covariance-PCA convention).
 *    2. y = PC1 score = (U * S).col(0) = X * v1  (ordinary PCA on centered X).
 *    3. lambda.min / lambda.1se from ridge_cv_glmnet(X, y)   [glmnet SD scale].
 *    4. lambda_2_lars = lambda.1se * p / 2                   [R's lm_dummy rule].
 *
 *  Why feeding centered-only X (instead of the L2-normalised X the selector
 *  uses internally) is correct here: ridge_cv_glmnet re-normalises every column
 *  to unit sample-SD inside each fold, which *overwrites* any prior per-column
 *  scaling. The selected lambda is therefore identical whether the input columns
 *  are raw, centered, or L2-normalised — and identical to R's cv.glmnet on the
 *  same centered X with standardize = TRUE. This lets us dump one shared
 *  (X, y) pair and compare all three (C++ glmnet path, C++ direct path, R).
 *
 *  For contrast the probe also runs ridge_cv_direct (the deactivated unit-L2
 *  path). Its lambda is expected to be ~(n - 1)x smaller than the glmnet path;
 *  the printed ratio makes the scale relationship explicit.
 *
 *  Outputs (written to DEMO_OUTPUT_DIR):
 *    - console + lambda2_probe_<snr>.txt : per-method lambda table
 *    - lambda2_probe_X.csv               : centered X for the dump case
 *    - lambda2_probe_y.csv               : PC1 score y for the dump case
 *  Cross-check in R with: R/.../trex_spca/lambda2_probe.R
 */
// ==============================================================================

// std includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/SVD>

// Ridge lambda-selection (header-only)
#include <ml_methods/model_selection/ridge_cv_glmnet.hpp>
#include <ml_methods/model_selection/ridge_cv_direct.hpp>
#include <ml_methods/model_selection/elastic_net_gaussian_strong_rule.hpp>

// Demo utilities (DataGenerator, FactorModelData, DEMO_OUTPUT_DIR consumers)
#include "trex_spca_sim_utils.hpp"

// ==============================================================================

using namespace spca_sim;
namespace ms = trex::ml_methods::model_selection;

// ==============================================================================
// Helpers
// ==============================================================================

/** @brief PC1 score y = (U * S).col(0) = X * v1 for a centered matrix X. */
static Eigen::VectorXd pc1_score(const Eigen::MatrixXd& X_centered) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        X_centered, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return svd.matrixU().col(0) * svd.singularValues()(0);
}

/** @brief Write a dense matrix to CSV (no header), full precision. */
static void dump_csv(const std::string& path, const Eigen::MatrixXd& M) {
    std::ofstream f(path);
    f << std::setprecision(17);
    for (Eigen::Index i = 0; i < M.rows(); ++i) {
        for (Eigen::Index j = 0; j < M.cols(); ++j) {
            f << M(i, j);
            if (j + 1 < M.cols()) f << ",";
        }
        f << "\n";
    }
}

/** @brief Per-method lambda summary for a single (X, y). */
struct LambdaProbe {
    double glmnet_min   = 0.0;
    double glmnet_1se   = 0.0;
    double glmnet_lars  = 0.0;   ///< lambda.1se * p / 2  (production lambda_2)
    double direct_min   = 0.0;
    double direct_1se   = 0.0;
    double direct_lars  = 0.0;
    double en_min       = 0.0;   ///< coordinate-descent ridge (alpha=0) lambda.min
    double en_1se       = 0.0;   ///< coordinate-descent ridge (alpha=0) lambda.1se
    double en_lars      = 0.0;   ///< en_1se * p / 2  (CD production lambda_2)
};

static LambdaProbe probe_lambda2(const Eigen::MatrixXd& X_centered,
                                 const Eigen::VectorXd& y,
                                 int          n_folds  = 10,
                                 Eigen::Index n_lambda = 100,
                                 unsigned     cv_seed  = 0)
{
    const double p = static_cast<double>(X_centered.cols());

    // glmnet-compatible (sample-SD normalisation, alpha = 0 lambda_max) ----
    ms::ridge_cv_glmnet cv_g;
    cv_g.fit(X_centered, y, n_folds, n_lambda, /*lambda_min_ratio=*/-1.0, cv_seed);

    // direct (unit-L2 normalisation, deactivated production path) -----------
    ms::ridge_cv_direct cv_d;
    cv_d.fit(X_centered, y, n_folds, n_lambda, /*lambda_ratio=*/1000.0, cv_seed);

    // coordinate-descent ridge (alpha=0), glmnet-faithful with fdev/devmax --
    // intercept=false to match the R cross-check (X and y are pre-centered) and
    // the production SPCA path (base-class centers X / y before computeLambda2).
    ms::elastic_net_cv_gaussian cv_en;
    cv_en.fit(X_centered, y, /*alpha=*/0.0, n_folds, n_lambda,
              /*lambda_min_ratio=*/-1.0, cv_seed,
              /*standardize=*/true, /*intercept=*/false);

    LambdaProbe out;
    out.glmnet_min  = cv_g.cv_min();
    out.glmnet_1se  = cv_g.cv_1se();
    out.glmnet_lars = out.glmnet_1se * p / 2.0;
    out.direct_min  = cv_d.cv_min();
    out.direct_1se  = cv_d.cv_1se();
    out.direct_lars = out.direct_1se * p / 2.0;
    out.en_min      = cv_en.cv_min();
    out.en_1se      = cv_en.cv_1se();
    out.en_lars     = out.en_1se * p / 2.0;
    return out;
}

// ==============================================================================
// Probe driver
// ==============================================================================

void run_lambda2_probe(const SimConfig&           cfg,
                       const std::vector<double>& snr_values,
                       std::size_t                n_seeds,
                       double                     dump_snr,
                       const std::string&         stem)
{
    const std::string folder = DEMO_OUTPUT_DIR;
    std::filesystem::create_directories(folder);
    std::ofstream out_file(folder + stem + ".txt");

    auto dual = [&](const std::string& s) {
        std::cout << s;
        if (out_file.is_open()) out_file << s;
    };

    {
        std::ostringstream ss;
        ss << "\n==============================================================\n"
           << "  lambda_2 determination probe  (PC1, n=" << cfg.n
           << ", p=" << cfg.p << ")\n"
           << "  averaged over " << n_seeds << " seeds per SNR\n"
           << "==============================================================\n\n"
           << std::left << std::setw(10) << "SNR(dB)"
           << std::right
           << std::setw(14) << "glmnet.1se"
           << std::setw(14) << "glmnet_lars"
           << std::setw(14) << "en.1se"
           << std::setw(14) << "en_lars"
           << std::setw(14) << "direct.1se"
           << std::setw(14) << "direct_lars"
           << std::setw(12) << "ratio"   // glmnet_lars / direct_lars ~ (n-1)
           << "\n"
           << std::string(106, '-') << "\n";
        dual(ss.str());
    }

    const double n_minus_1 = static_cast<double>(cfg.n - 1);

    for (double snr : snr_values) {
        double acc_g1se = 0.0, acc_glars = 0.0;
        double acc_d1se = 0.0, acc_dlars = 0.0, acc_ratio = 0.0;
        double acc_en1se = 0.0, acc_enlars = 0.0;

        for (std::size_t s = 0; s < n_seeds; ++s) {
            const unsigned seed =
                static_cast<unsigned>(cfg.base_seed) + static_cast<unsigned>(s) * 1000u
                + static_cast<unsigned>(static_cast<int>(std::round(snr)));

            auto dat = DataGenerator::generate_sparse_factor_model(
                cfg.n, cfg.p, cfg.p1, cfg.M, snr, cfg.overlap_pool_size, seed);

            // Pipeline preprocessing: center X once (covariance-PCA convention).
            dat.X.rowwise() -= dat.X.colwise().mean();
            const Eigen::VectorXd y = pc1_score(dat.X);

            const LambdaProbe lp = probe_lambda2(dat.X, y);
            acc_g1se  += lp.glmnet_1se;
            acc_glars += lp.glmnet_lars;
            acc_d1se  += lp.direct_1se;
            acc_dlars += lp.direct_lars;
            acc_en1se  += lp.en_1se;
            acc_enlars += lp.en_lars;
            acc_ratio += (lp.direct_lars > 0.0)
                       ? lp.glmnet_lars / lp.direct_lars : 0.0;

            // Dump one representative (snr, first seed) for the R cross-check.
            if (s == 0 && std::round(snr) == std::round(dump_snr)) {
                dump_csv(folder + stem + "_X.csv", dat.X);
                dump_csv(folder + stem + "_y.csv", y);
                std::ostringstream ds;
                ds << "  [dump] SNR=" << snr << " seed=" << seed
                   << " -> " << stem << "_X.csv / " << stem << "_y.csv\n"
                   << "  [dump single-case lambdas]"
                   << "  glmnet.min=" << std::setprecision(5) << lp.glmnet_min
                   << "  glmnet.1se=" << lp.glmnet_1se
                   << "  en.min=" << lp.en_min
                   << "  en.1se=" << lp.en_1se << "\n"
                   << "                            "
                   << "  (compare to R: lambda.min=97.64, lambda.1se=151.5+/-17)\n";
                dual(ds.str());
            }
        }

        const double inv = 1.0 / static_cast<double>(n_seeds);
        std::ostringstream row;
        row << std::left << std::setw(10) << std::fixed << std::setprecision(1) << snr
            << std::right << std::setprecision(5)
            << std::setw(14) << acc_g1se  * inv
            << std::setw(14) << acc_glars * inv
            << std::setw(14) << acc_en1se  * inv
            << std::setw(14) << acc_enlars * inv
            << std::setw(14) << acc_d1se  * inv
            << std::setw(14) << acc_dlars * inv
            << std::setprecision(2)
            << std::setw(12) << acc_ratio * inv
            << "\n";
        dual(row.str());
    }

    {
        std::ostringstream ss;
        ss << std::string(106, '-') << "\n"
           << "  Expected ratio glmnet_lars / direct_lars  ~  n - 1 = "
           << std::fixed << std::setprecision(1) << n_minus_1 << "\n"
           << "  production lambda_2 = glmnet_lars column (CV_1SE_GLMNET * p/2)\n"
           << "  en_lars = coordinate-descent ridge (CV_1SE_EN * p/2); should\n"
           << "  track glmnet_lars (CD lambda.min matches cv.glmnet exactly).\n\n";
        dual(ss.str());
    }

    if (out_file.is_open()) {
        std::cout << "[Info] Probe table saved to: " << folder + stem + ".txt\n";
        std::cout << "[Info] Cross-check in R:      lambda2_probe.R "
                  << "(loads " << stem << "_X.csv / " << stem << "_y.csv)\n\n";
        out_file.close();
    }
}

// ==============================================================================
// Main
// ==============================================================================

int main() {
    SimConfig cfg;          // n=50, p=100, p1=5, M=3, overlap_pool_size=30, base_seed=42
    cfg.num_MC = 0;         // unused by the probe

    const std::vector<double> snr_values = {-10.0, -7.0, -5.0, -3.0, 0.0, 5.0, 10.0};

    run_lambda2_probe(cfg, snr_values,
                      /*n_seeds=*/20,
                      /*dump_snr=*/-10.0,   // dump the problematic low-SNR case for R
                      /*stem=*/"lambda2_probe");

    return 0;
}
