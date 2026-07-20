// ==============================================================================
// demo_table_utils.hpp
// ==============================================================================
/**
 * @file demo_table_utils.hpp
 *
 * @brief Shared layout constants and column chunking for the demo result tables.
 *
 * @details
 *  Every suite prints its Monte Carlo summaries in the same layout: the series
 *  name on its own line, metric rows indented beneath it, and the sweep axis as
 *  a single header row.
 *
 *      SNR                          0.10      0.20      0.30
 *    --------------------------------------------------------
 *
 *    TLARS
 *      FDR                        0.0200    0.0363    0.0437
 *      TPR                        0.0050    0.2830    0.7395
 *
 *  A wide sweep grid would push that past the terminal width and wrap at an
 *  arbitrary point, which destroys the column alignment. column_chunks() splits
 *  the sweep into consecutive blocks that each fit kMaxLineWidth; callers emit
 *  one header + one set of metric rows per block, repeating the series name
 *  with a "(continued)" marker on the second and later blocks.
 *
 *  Demos with a dense grid (e.g. the 21-point SNR sweeps of trex demos 03 and
 *  04) therefore stay readable without thinning the grid, which the figures
 *  need at full resolution.
 *
 * @note Demo-internal header — not part of the TRexSelector library.
 */
// ==============================================================================

#ifndef TREX_DEMOS_DEMO_TABLE_UTILS_HPP
#define TREX_DEMOS_DEMO_TABLE_UTILS_HPP

// std includes
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace demo_tables {

/** @brief Leading indent of a metric row. */
inline constexpr int kIndentW = 2;

/** @brief Width of the left-aligned metric-label column. */
inline constexpr int kMetricW = 23;

/** @brief Longest line a table may produce before the sweep is split. */
inline constexpr std::size_t kMaxLineWidth = 120;

/**
 * @brief Split a sweep axis into consecutive column blocks that each fit.
 *
 * @param n_cols  Number of sweep values.
 * @param col_w   Width of one value column.
 * @return Half-open [begin, end) index ranges, one per block. Always at least
 *         one block, even if a single column already exceeds the budget.
 */
inline std::vector<std::pair<std::size_t, std::size_t>>
column_chunks(std::size_t n_cols, int col_w)
{
    const std::size_t label_w =
        static_cast<std::size_t>(kIndentW + kMetricW);
    std::size_t per_chunk = 1;
    if (col_w > 0 && kMaxLineWidth > label_w) {
        per_chunk = (kMaxLineWidth - label_w) / static_cast<std::size_t>(col_w);
        if (per_chunk == 0) per_chunk = 1;
    }

    std::vector<std::pair<std::size_t, std::size_t>> chunks;
    for (std::size_t b = 0; b < n_cols; b += per_chunk) {
        chunks.emplace_back(b, std::min(b + per_chunk, n_cols));
    }
    if (chunks.empty()) chunks.emplace_back(0u, 0u);
    return chunks;
}

/**
 * @brief Series heading for a given block: plain name, then "(continued)".
 *
 * @param name        Series (solver / method) name.
 * @param chunk_index Zero-based block index.
 */
inline std::string series_heading(const std::string& name,
                                  std::size_t chunk_index)
{
    return (chunk_index == 0) ? name : name + "  (continued)";
}

} // namespace demo_tables

// ==============================================================================
#endif /* End of TREX_DEMOS_DEMO_TABLE_UTILS_HPP */
// ==============================================================================
