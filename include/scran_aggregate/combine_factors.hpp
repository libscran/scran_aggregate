#ifndef SCRAN_AGGREGATE_COMBINE_FACTORS_HPP
#define SCRAN_AGGREGATE_COMBINE_FACTORS_HPP

#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>

/**
 * @file combine_factors.hpp
 * @brief Combine categorical factors into a single factor. 
 */

namespace scran_aggregate {

/**
 * @tparam Factor_ Factor type.
 * Any type may be used here as long as it is comparable.
 * @tparam Combined_ Integer type for the combined factor.
 * This should be large enough to hold the number of unique combinations.
 *
 * @param n Number of observations (i.e., cells).
 * @param[in] factors Vector of pointers to arrays of length `n`, each containing a different factor.
 * @param[out] combined Pointer to an array of length `n` in which the combined factor is to be stored.
 * On output, each entry determines the corresponding observation's combination of levels by indexing into the inner vectors of the returned object,
 * i.e., `j := combined[i]` represents the combination `(output[0][j], output[1][j], ...)`.
 *
 * @return 
 * Vector of vectors containing each unique combinations of factor levels.
 * Each inner vector corresponds to a factor in `factors`, and all inner vectors have the same length.
 * Corresponding entries of the inner vectors define a particular combination of levels,
 * i.e., the first combination is defined as `(output[0][0], output[1][0], ...)`,
 * the second combination is defined as `(output[0][1], output[1][1], ...)`, and so on.
 * Combinations are guaranteed to be sorted by the first factor, then the second, etc.
 */
template<typename Factor_, typename Combined_>
std::vector<std::vector<Factor_> > combine_factors(size_t n, const std::vector<const Factor_*>& factors, Combined_* combined) {
    // Using a map with a custom comparator that uses the index
    // of first occurrence of each factor as the key. Currently using a map
    // to (i) avoid issues with collisions of combined hashes and (ii)
    // avoid having to write more code for sorting a vector of arrays.
    auto cmp = [&](size_t left, size_t right) -> bool {
        for (auto curf : factors) {
            if (curf[left] < curf[right]) {
                return true;
            } else if (curf[left] > curf[right]) {
                return false;
            }
        }
        return false;
    };

    auto eq = [&](size_t left, size_t right) -> bool {
        for (auto curf : factors) {
            if (curf[left] != curf[right]) {
                return false;
            }
        }
        return true;
    };

    std::map<size_t, Combined_, decltype(cmp)> mapping(cmp);
    for (size_t i = 0; i < n; ++i) {
        auto mIt = mapping.find(i);
        if (mIt == mapping.end() || !eq(i, mIt->first)) {
            mapping.insert(mIt, std::pair<size_t, Combined_>(i, 0));
        }
    }

    // Obtaining the sorted set of unique combinations; easy to do for a
    // map because it's already sorted!
    size_t nfac = factors.size();
    std::vector<std::vector<Factor_> > output(nfac);
    size_t nuniq = mapping.size();
    for (auto& ofac : output) {
        ofac.reserve(nuniq);
    }

    auto mIt = mapping.begin();
    for (size_t u = 0; u < nuniq; ++u, ++mIt) {
        auto ix = mIt->first;
        for (size_t f = 0; f < nfac; ++f) {
            output[f].push_back(factors[f][ix]);
        }
        mIt->second = u;
    }

    // Mapping each cell to its unique combination.
    for (size_t i = 0; i < n; ++i) {
        combined[i] = mapping[i];
    }

    return output;
}

}

#endif
