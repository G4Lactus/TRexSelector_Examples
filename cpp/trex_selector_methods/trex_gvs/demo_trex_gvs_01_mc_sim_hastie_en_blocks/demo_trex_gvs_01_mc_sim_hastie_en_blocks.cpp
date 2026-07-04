// ==============================================================================
// demo_trex_gvs_01_mc_sim_hastie_en_blocks.cpp
// ==============================================================================
/**
 * @file demo_trex_gvs_01_mc_sim_hastie_en_blocks.cpp
 *
 * @brief MC simulation: T-Rex+GVS on the Hastie DGP.
 *  This scenario evaluates grouped-signal recovery under strong equicorrelation,
 *  not sparse within-group support identification.
 *  Within-group correlation is controlled through sd_x, with rho = 1 / (1 + sd_x^2).
 *
 *
 * @details
 *  Five parts, each in its own function:
 *    sim_hastie_part1() — SNR sweep (fixed sd_x = sqrt(0.01))
 *    sim_hastie_part2() — rho sweep (fixed SNR = 2.0, sd_x derived from rho)
 *    sim_hastie_part3() — 2-D SNR × rho grid
 *    sim_hastie_part4() — lambda2_method comparison: CV_1SE_SVD / CV_1SE_CCD
 *    sim_hastie_part5() — hc_linkage comparison: Single / Complete / Average / WPGMA
 *
 * DGP: three equicorrelated groups of 50 variables each, all active (s = 150),
 *      plus 350 null white-noise variables; p = 500 and n = 200.
 *  Selector: K=20, tFDR=0.1, corr_max=0.98, lambda2=CV_1SE_CCD.
 *  Methods compared in Parts 1–3: EN (TENET), EN (TENET_AUG), IEN.
 *  MC: 200 trials per grid point.
 */
// ==============================================================================

#include "trex_gvs_dgps.hpp"
#include "trex_gvs_simulation_utils.hpp"


// ==============================================================================

// Namespace usage
using namespace gvs_demo;


// ==============================================================================
//  Part 1 — SNR sweep
// ==============================================================================

/**
 * @brief Run an MC SNR sweep on the Hastie DGP at fixed within-group correlation.
 *
 * Uses sd_x = sqrt(0.01) (rho ≈ 0.99) and compares EN (TENET), EN (TENET_AUG),
 * and IEN across the SNR grid; prints and saves mean FDP/TPP tables.
 */

static void sim_hastie_part1(const GVSSimConfig& cfg, const std::string& tag)
{
    const std::vector<double> snr_grid = {
        0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    GVSDGPFactory snr_fn = [&cfg](double snr, unsigned seed) {
        return make_hastie_dgp(cfg.n, cfg.p, snr, cfg.sd_x, seed);
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
        "  s=150  MC=" + std::to_string(cfg.num_MC) +
        "  tFDR=" + fmt_num(cfg.tFDR) +
        "  sd_x=sqrt(0.01)");

    auto en_snr     = run_gvs_snr_sweep(
        snr_fn, snr_grid, cfg, gvs_ctrl_en,     trex_ctrl, "EN");

    auto en_aug_snr = run_gvs_snr_sweep(
        snr_fn, snr_grid, cfg, gvs_ctrl_en_aug, trex_ctrl, "EN+AUG");

    auto ien_snr    = run_gvs_snr_sweep(
        snr_fn, snr_grid, cfg, gvs_ctrl_ien,    trex_ctrl, "IEN");

    print_mc_snr_table(tag, snr_grid, en_snr, en_aug_snr, ien_snr,
                       cfg, "gvs_" + tag + "_snr");
}



// ==============================================================================
//  Part 2 — rho sweep
// ==============================================================================

/**
 * @brief Run an MC rho sweep on the Hastie DGP at fixed SNR.
 *
 * Converts each target rho to sd_x = sqrt((1 - rho) / rho), then compares
 * EN (TENET), EN (TENET_AUG), and IEN across the rho grid.
 */

static void sim_hastie_part2(const GVSSimConfig& cfg, const std::string& tag)
{
    const std::vector<double> rho_grid = {
        0.10, 0.20, 0.30, 0.50, 0.70,
        0.80, 0.90, 0.95, 0.99};

    GVSRhoDGPFactory rho_fn = [&cfg](double rho, unsigned seed) {
        const double sd_x = std::sqrt((1.0 - rho) / rho);
        return make_hastie_dgp(cfg.n, cfg.p, cfg.snr, sd_x, seed);
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
        "  (SNR=" + fmt_num(cfg.snr) + ")");

    auto en_rho     = run_gvs_rho_sweep(
        rho_fn, rho_grid, cfg, gvs_ctrl_en,     trex_ctrl, "EN");

    auto en_aug_rho = run_gvs_rho_sweep(
        rho_fn, rho_grid, cfg, gvs_ctrl_en_aug, trex_ctrl, "EN+AUG");

    auto ien_rho    = run_gvs_rho_sweep(
        rho_fn, rho_grid, cfg, gvs_ctrl_ien,    trex_ctrl, "IEN");

    print_mc_rho_table(tag, rho_grid, en_rho, en_aug_rho,
                       ien_rho, cfg, "gvs_" + tag + "_rho");
}



// ==============================================================================
//  Part 3 — 2-D SNR × rho sweep
// ==============================================================================

/**
 * @brief Run a 2-D MC sweep over SNR and rho on the Hastie DGP.
 *
 * Evaluates EN (TENET), EN (TENET_AUG), and IEN on an SNR × rho grid, then
 * prints and saves matrix summaries for mean TPP and mean FDP.
 */

static void sim_hastie_part3(const GVSSimConfig& cfg, const std::string& tag)
{
    const std::vector<double> snr_grid_2d = {
        0.2, 0.5, 1.0, 2.0, 5.0};
    const std::vector<double> rho_grid_2d = {
        0.30, 0.50, 0.70, 0.90, 0.95, 0.99};

    GVS2DDGPFactory dgp_2d = [&cfg](double snr, double rho, unsigned seed) {
        const double sd_x = std::sqrt((1.0 - rho) / rho);
        return make_hastie_dgp(cfg.n, cfg.p, snr, sd_x, seed);
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

    auto en_2d     = run_gvs_2d_sweep(
        dgp_2d, snr_grid_2d, rho_grid_2d, cfg,
        gvs_ctrl_en,     trex_ctrl, "EN");

    auto en_aug_2d = run_gvs_2d_sweep(
        dgp_2d, snr_grid_2d, rho_grid_2d, cfg,
        gvs_ctrl_en_aug, trex_ctrl, "EN+AUG");

    auto ien_2d    = run_gvs_2d_sweep(
        dgp_2d, snr_grid_2d, rho_grid_2d, cfg,
        gvs_ctrl_ien,    trex_ctrl, "IEN");

    std::vector<std::string> snr_labels, rho_labels;
    for (double s : snr_grid_2d)
        snr_labels.push_back("snr=" + fmt_num(s));
    for (double r : rho_grid_2d)
        rho_labels.push_back("rho=" + fmt_num(r));

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
    save_mc_2d_tables(tag, snr_labels, rho_labels, en_2d,
                      en_aug_2d, ien_2d, cfg, "gvs_" + tag + "_2d");
}



// ==============================================================================
//  Part 4 — lambda2_method comparison
// ==============================================================================

/**
 * @brief Compare lambda2 selection methods at a fixed Hastie operating point.
 *
 * Runs EN and IEN at SNR = 2.0 and rho = 0.7, comparing CV_1SE_SVD against
 * CV_1SE_CCD while holding the rest of the selector settings fixed.
 */

static void sim_hastie_part4(const GVSSimConfig& cfg, const std::string& tag)
{
    // Fixed operating point: SNR=2.0, rho=0.7.
    // sd_x = sqrt((1 - rho) / rho) = sqrt(0.3 / 0.7)
    constexpr double fixed_snr = 2.0;
    constexpr double fixed_rho = 0.7;
    const double     fixed_sdx = std::sqrt((1.0 - fixed_rho) / fixed_rho);

    GVSDGPFactory dgp_fn = [&cfg, fixed_snr, fixed_sdx](double /*unused*/, unsigned seed) {
        return make_hastie_dgp(cfg.n, cfg.p, fixed_snr, fixed_sdx, seed);
    };

    auto trex_ctrl = make_gvs_trex_control(cfg.K);

    const std::vector<LambdaSelectionMethod> methods = {
        LambdaSelectionMethod::CV_1SE_SVD,
        LambdaSelectionMethod::CV_1SE_CCD
    };
    const std::vector<std::string> method_labels = {"CV_1SE_SVD", "CV_1SE_CCD"};

    cdiag::print_section_header(
        "Part 4: lambda2_method Comparison  |  " + tag +
        "\nn=" + std::to_string(cfg.n) + "  p=" + std::to_string(cfg.p) +
        "  SNR=" + fmt_num(fixed_snr) +
        "  rho=" + fmt_num(fixed_rho) +
        "  s=150  MC=" + std::to_string(cfg.num_MC) +
        "  tFDR=" + fmt_num(cfg.tFDR) +
        "  methods: CV_1SE_SVD / CV_1SE_CCD");

    std::vector<GVSGridPointResult> en_results, ien_results;
    en_results.reserve(methods.size());
    ien_results.reserve(methods.size());

    for (std::size_t i = 0; i < methods.size(); ++i) {
        TRexGVSControlParameter gvs_ctrl_en;
        gvs_ctrl_en.gvs_type       = GVSType::EN;
        gvs_ctrl_en.lambda2_method = methods[i];
        gvs_ctrl_en.corr_max       = cfg.corr_max;
        gvs_ctrl_en.hc_linkage     = hac::LinkageMethod::Single;
        gvs_ctrl_en.en_solver      = ENSolverType::TENET;

        en_results.push_back(
            run_gvs_mc_trials(dgp_fn, 0.0, cfg.num_MC,
                              "EN+" + method_labels[i],
                              cfg.tFDR, gvs_ctrl_en, trex_ctrl,
                              static_cast<unsigned>(cfg.base_seed)));

        TRexGVSControlParameter gvs_ctrl_ien;
        gvs_ctrl_ien.gvs_type       = GVSType::IEN;
        gvs_ctrl_ien.lambda2_method = methods[i];
        gvs_ctrl_ien.corr_max       = cfg.corr_max;
        gvs_ctrl_ien.hc_linkage     = hac::LinkageMethod::Single;

        ien_results.push_back(
            run_gvs_mc_trials(dgp_fn, 0.0, cfg.num_MC,
                                 "IEN+" + method_labels[i],
                                 cfg.tFDR,
                                 gvs_ctrl_ien, trex_ctrl,
                                 static_cast<unsigned>(cfg.base_seed)));
    }

    {
        const std::string dir(DEMO_OUTPUT_DIR);
        const std::string stem = "gvs_" + tag + "_lambda2_method";
        std::filesystem::create_directories(dir);
        std::ofstream txt_file(dir + stem + ".txt");
        std::ofstream* fout = txt_file.is_open() ? &txt_file : nullptr;

        auto write = [&](const std::string& s) {
            std::cout << s;
            if (fout) *fout << s;
        };

        {
            std::ostringstream hdr;
            hdr << "\n" << std::string(78, '=') << "\n"
                << "  T-Rex+GVS MC: " << tag << " — lambda2_method comparison"
                << "  (SNR=" << fixed_snr << ", rho=" << fixed_rho
                << ", MC=" << cfg.num_MC << ", tFDR=" << cfg.tFDR << ")\n"
                << std::string(78, '=') << "\n";
            write(hdr.str());
        }

        print_mc_1d_method_table("EN",  method_labels, en_results,  fout);
        print_mc_1d_method_table("IEN", method_labels, ien_results, fout);

        write(std::string(78, '=') + "\n\n");

        if (txt_file.is_open()) {
            txt_file.close();
            std::cout << "[Info] TXT saved to: " << dir + stem + ".txt\n";
        }
        std::ofstream csv(dir + stem + ".csv");
        if (csv.is_open()) {
            csv << "method,metric,lambda2_method,mean,sd\n"
                << std::fixed << std::setprecision(6);
            for (std::size_t i = 0; i < method_labels.size(); ++i) {
                const std::string& lbl = method_labels[i];
                csv << "EN,mean_FDP,"  << lbl << "," << en_results[i].mean_fdp  << "," << en_results[i].sd_fdp  << "\n";
                csv << "EN,mean_TPP,"  << lbl << "," << en_results[i].mean_tpp  << "," << en_results[i].sd_tpp  << "\n";
                csv << "IEN,mean_FDP," << lbl << "," << ien_results[i].mean_fdp << "," << ien_results[i].sd_fdp << "\n";
                csv << "IEN,mean_TPP," << lbl << "," << ien_results[i].mean_tpp << "," << ien_results[i].sd_tpp << "\n";
            }
            std::cout << "[Info] CSV saved to: " << dir + stem + ".csv\n\n";
        }
    }
}


// ==============================================================================
//  Part 5 — hc_linkage comparison
// ==============================================================================

/**
 * @brief Compare HAC linkage rules at a fixed Hastie operating point.
 *
 * Runs EN and IEN at SNR = 2.0 and rho = 0.7, comparing Single, Complete,
 * Average, and WPGMA linkage with lambda2 fixed at CV_1SE_CCD.
 */

static void sim_hastie_part5(const GVSSimConfig& cfg, const std::string& tag)
{
    // Same fixed operating point as Part 4: SNR=2.0, rho=0.7.
    // CV_1SE_CCD used throughout to isolate the linkage effect from lambda selection.
    constexpr double fixed_snr = 2.0;
    constexpr double fixed_rho = 0.7;
    const double     fixed_sdx = std::sqrt((1.0 - fixed_rho) / fixed_rho);

    GVSDGPFactory dgp_fn = [&cfg, fixed_snr, fixed_sdx](double /*unused*/, unsigned seed) {
        return make_hastie_dgp(cfg.n, cfg.p, fixed_snr, fixed_sdx, seed);
    };

    auto trex_ctrl = make_gvs_trex_control(cfg.K);

    const std::vector<hac::LinkageMethod> linkages = {
        hac::LinkageMethod::Single,
        hac::LinkageMethod::Complete,
        hac::LinkageMethod::Average,
        hac::LinkageMethod::WPGMA
    };
    const std::vector<std::string> linkage_labels = {
        "Single", "Complete", "Average", "WPGMA"};

    cdiag::print_section_header(
        "Part 5: hc_linkage Comparison  |  " + tag +
        "\nn=" + std::to_string(cfg.n) + "  p=" + std::to_string(cfg.p) +
        "  SNR=" + fmt_num(fixed_snr) +
        "  rho=" + fmt_num(fixed_rho) +
        "  s=150  MC=" + std::to_string(cfg.num_MC) +
        "  tFDR=" + fmt_num(cfg.tFDR) +
        "  lambda2=CV_1SE_CCD");

    std::vector<GVSGridPointResult> en_results, ien_results;
    en_results.reserve(linkages.size());
    ien_results.reserve(linkages.size());

    for (std::size_t i = 0; i < linkages.size(); ++i) {
        TRexGVSControlParameter gvs_ctrl_en;
        gvs_ctrl_en.gvs_type       = GVSType::EN;
        gvs_ctrl_en.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
        gvs_ctrl_en.corr_max       = cfg.corr_max;
        gvs_ctrl_en.hc_linkage     = linkages[i];
        gvs_ctrl_en.en_solver      = ENSolverType::TENET;

        en_results.push_back(
            run_gvs_mc_trials(dgp_fn, 0.0, cfg.num_MC,
                              "EN+" + linkage_labels[i],
                              cfg.tFDR, gvs_ctrl_en, trex_ctrl,
                              static_cast<unsigned>(cfg.base_seed)));

        TRexGVSControlParameter gvs_ctrl_ien;
        gvs_ctrl_ien.gvs_type       = GVSType::IEN;
        gvs_ctrl_ien.lambda2_method = LambdaSelectionMethod::CV_1SE_CCD;
        gvs_ctrl_ien.corr_max       = cfg.corr_max;
        gvs_ctrl_ien.hc_linkage     = linkages[i];

        ien_results.push_back(
            run_gvs_mc_trials(dgp_fn, 0.0, cfg.num_MC,
                              "IEN+" + linkage_labels[i],
                              cfg.tFDR, gvs_ctrl_ien, trex_ctrl,
                              static_cast<unsigned>(cfg.base_seed)));
    }

    {
        const std::string dir(DEMO_OUTPUT_DIR);
        const std::string stem = "gvs_" + tag + "_hc_linkage";
        std::filesystem::create_directories(dir);
        std::ofstream txt_file(dir + stem + ".txt");
        std::ofstream* fout = txt_file.is_open() ? &txt_file : nullptr;

        auto write = [&](const std::string& s) {
            std::cout << s;
            if (fout) *fout << s;
        };

        {
            std::ostringstream hdr;
            hdr << "\n" << std::string(78, '=') << "\n"
                << "  T-Rex+GVS MC: " << tag << " — hc_linkage comparison"
                << "  (SNR=" << fixed_snr << ", rho=" << fixed_rho
                << ", MC=" << cfg.num_MC << ", tFDR=" << cfg.tFDR << ")\n"
                << std::string(78, '=') << "\n";
            write(hdr.str());
        }

        print_mc_1d_method_table("EN",  linkage_labels, en_results,  fout);
        print_mc_1d_method_table("IEN", linkage_labels, ien_results, fout);

        write(std::string(78, '=') + "\n\n");

        if (txt_file.is_open()) {
            txt_file.close();
            std::cout << "[Info] TXT saved to: " << dir + stem + ".txt\n";
        }
        std::ofstream csv(dir + stem + ".csv");
        if (csv.is_open()) {
            csv << "method,metric,hc_linkage,mean,sd\n"
                << std::fixed << std::setprecision(6);
            for (std::size_t i = 0; i < linkage_labels.size(); ++i) {
                const std::string& lbl = linkage_labels[i];
                csv << "EN,mean_FDP,"  << lbl << "," << en_results[i].mean_fdp  << "," << en_results[i].sd_fdp  << "\n";
                csv << "EN,mean_TPP,"  << lbl << "," << en_results[i].mean_tpp  << "," << en_results[i].sd_tpp  << "\n";
                csv << "IEN,mean_FDP," << lbl << "," << ien_results[i].mean_fdp << "," << ien_results[i].sd_fdp << "\n";
                csv << "IEN,mean_TPP," << lbl << "," << ien_results[i].mean_tpp << "," << ien_results[i].sd_tpp << "\n";
            }
            std::cout << "[Info] CSV saved to: " << dir + stem + ".csv\n\n";
        }
    }
}


// ==============================================================================
//  main
// ==============================================================================

int main()
{

    // ==========================================================================
    //  Simulation parameters
    // ==========================================================================

    std::cout.setf(std::ios::unitbuf);
    omp_set_num_threads(6);
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    GVSSimConfig cfg;
    cfg.n         = 200;
    cfg.p         = 500;
    cfg.sd_x      = 0.1;    // sqrt(0.01) => rho ≈ 0.99
    cfg.tFDR      = 0.1;
    cfg.K         = 20;
    cfg.num_MC    = 200;
    cfg.base_seed = 2026;
    cfg.corr_max  = 0.98;
    cfg.snr       = 2.0;    // fixed SNR for Part 2 (rho sweep)

    const std::string scenario_tag = "Hastie";

    // ==========================================================================
    //  Run simulations
    // ==========================================================================

    if (true) sim_hastie_part1(cfg, scenario_tag);
    if (true) sim_hastie_part2(cfg, scenario_tag);
    if (true) sim_hastie_part3(cfg, scenario_tag);
    if (true) sim_hastie_part4(cfg, scenario_tag);
    if (true) sim_hastie_part5(cfg, scenario_tag);

    std::cout << scenario_tag << " GVS MC simulations complete.\n";
    return 0;
}
