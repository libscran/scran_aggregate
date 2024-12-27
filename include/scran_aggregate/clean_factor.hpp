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
 * Any type may be used here as long as it is hashable and has an equality operator.
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
    auto unique = [&]{ // scoping this in an IIFE to release map memory sooner.
        std::unordered_map<Factor_, Output_> mapping;
        for (size_t i = 0; i < n; ++i) {
            auto current = factor[i];
            auto mIt = mapping.find(current);
            if (mIt != mapping.end()) {
                cleaned[i] = mIt->second;
            } else {
                Output_ alt = mapping.size();
                mapping[current] = alt;
                cleaned[i] = alt;
            }
        }
        return std::vector<std::pair<Factor_, Output_> >(mapping.begin(), mapping.end());
    }();

    // Remapping to a sorted set.
    std::sort(unique.begin(), unique.end());
    size_t nuniq = unique.size();
    std::vector<Output_> remapping(nuniq);
    std::vector<Factor_> output(nuniq);
    for (size_t u = 0; u < nuniq; ++u) {
        remapping[unique[u].second] = u;
        output[u] = unique[u].first;
    }

    // Mapping each cell to its sorted factor.
    for (size_t i = 0; i < n; ++i) {
        cleaned[i] = remapping[cleaned[i]];
    }

    return output;
}

}

#endif
