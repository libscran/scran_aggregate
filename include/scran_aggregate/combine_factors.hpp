#ifndef SCRAN_AGGREGATE_COMBINE_FACTORS_HPP
#define SCRAN_AGGREGATE_COMBINE_FACTORS_HPP

#include <algorithm>
#include <vector>
#include <map>

/**
 * @file combine_factors.hpp
 * @brief Combine categorical factors into a single factor. 
 */

namespace scran_aggregate {

/**
 * @brief Unique combinations from `combine_factors()`.
 *
 * @tparam Factor_ Factor type, typically an integer.
 */
template<typename Factor_>
struct FactorCombinations {
    /**
     * @cond
     */
    FactorCombinations(size_t n) : factors(n) {}
    /**
     * @endcond
     */

    /**
     * Unique combinations of factor levels.
     * Each inner vector corresponds to a factor used as input to `combine_factors()`.
     * All inner vectors have the same length.
     * Corresponding entries of the inner vectors define a particular combination of levels,
     * i.e., the first combination is defined as `(factors[0][0], factors[1][0], ...)`,
     * the second combination is defined as `(factors[0][1], factors[1][1], ...)`, and so on.
     * Combinations are guaranteed to be sorted.
     */
    std::vector<std::vector<Factor_> > factors;

    /**
     * Number of cells in each unique combination of factor levels.
     * This has the same length as each inner vector of `factors`.
     * All entries are guaranteed to be positive.
     */
    std::vector<size_t> counts;
};

/**
 * @tparam Factor_ Factor type.
 * Any type may be used here as long as it is comparable.
 * @tparam Combined_ Integer type for the combined factor.
 * This should be large enough to hold the number of unique combinations.
 *
 * @param n Number of observations (i.e., cells).
 * @param[in] factors Pointers to arrays of length `n`, each containing a different factor.
 * @param[out] combined Pointer to an array of length `n`, in which the combined factor is to be stored.
 *
 * @return 
 * Object containing the unique combinations of levels observed in `factors`.
 * The combined factor written to `combined`, where each entry is an index into the relevant combination of the output `FactorCombinations` object.
 */
template<typename Factor_, typename Combined_>
FactorCombinations<Factor_> combine_factors(size_t n, const std::vector<const Factor_*>& factors, Combined_* combined) {
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

    // Obtaining the sorted set of unique combinations.
    size_t nfac = factors.size();
    FactorCombinations<Factor_> output(nfac);
    size_t nuniq = mapping.size();
    for (auto& ofac : output.factors) {
        ofac.reserve(nuniq);
    }
    output.counts.resize(nuniq);

    auto mIt = mapping.begin();
    for (size_t u = 0; u < nuniq; ++u, ++mIt) {
        auto ix = mIt->first;
        for (size_t f = 0; f < nfac; ++f) {
            output.factors[f].push_back(factors[f][ix]);
        }
        mIt->second = u;
    }

    // Mapping each cell to its unique combination.
    for (size_t i = 0; i < n; ++i) {
        auto chosen = mapping[i];
        combined[i] = chosen;
        ++(output.counts[chosen]);
    }

    return output;
}

}

#endif
