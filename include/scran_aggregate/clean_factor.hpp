#ifndef SCRAN_AGGREGATE_CLEAN_FACTORS_HPP
#define SCRAN_AGGREGATE_CLEAN_FACTORS_HPP

#include <unordered_map>
#include <vector>
#include <algorithm>

/**
 * @file clean_factor.hpp
 * @brief Clean up a categorical factor.
 */

namespace scran_aggregate {

/**
 * Clean up a categorical factor by removing unused levels. 
 * This yields the same results as `combine_factors()` with a single factor.
 *
 * @tparam Factor_ Factor type.
 * Any type may be used here as long as it is comparable.
 * @tparam Output_ Integer type for the cleaned factor.
 *
 * @param n Number of observations (i.e., cells).
 * @param[in] factor Pointer to an array of length `n` containing a factor.
 * @param[out] cleaned Pointer to an array of length `n` in which the cleaned factor is to be stored.
 * All values are integers in \f$[0, N)\f$ where \f$N\f$ is the length of the output vector;
 * all integers in this range are guaranteed to be present at least once in `cleaned`.
 *
 * @return A sorted vector of the original levels that were observed at least once in `factor`.
 * For any observation `i`, it is guaranteed that `output[cleaned[i]] == factor[i]`.
 */
template<typename Factor_, typename Output_>
std::vector<Factor_> clean_factor(size_t n, const Factor_* factor, Output_* cleaned) {
    std::unordered_map<Factor_, Output_> mapping;
    for (size_t i = 0; i < n; ++i) {
        auto current = factor[i];
        mapping[current] = 0;
    }

    // Obtaining the sorted set of unique combinations.
    std::vector<Factor_> output;
    size_t nuniq = mapping.size();
    output.reserve(nuniq);
    for (const auto& mp : mapping) {
        output.push_back(mp.first);
    }
    std::sort(output.begin(), output.end());

    Output_ counter = 0;
    for (auto key : output) {
        mapping[key] = counter;
        ++counter;
    }

    // Mapping each cell to its sorted factor.
    for (size_t i = 0; i < n; ++i) {
        cleaned[i] = mapping[factor[i]];
    }

    return output;
}

}

#endif
