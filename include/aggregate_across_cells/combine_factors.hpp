#ifndef SCRAN_AGGREGATE_ACROSS_CELLS_COMBINE_FACTORS_HPP
#define SCRAN_AGGREGATE_ACROSS_CELLS_COMBINE_FACTORS_HPP

#include <algorithm>
#include <vector>
#include <map>

/**
 * @file combine_factors.hpp
 * @brief Combine categorical factors into a single factor. 
 */

namespace aggregate_across_cells {

/**
 * @brief Unique combinations of factors.
 *
 * @tparam Factor Factor type, typically an integer.
 */
template<typename Factor>
struct Combinations {
    /**
     * @cond
     */
    Combinations(size_t n) : factors(n) {}
    /**
     * @endcond
     */

    /**
     * Unique combinations of factor levels.
     * Each inner vector corresponds to a factor.
     * All inner vectors have the same length.
     * Corresponding entries of the inner vectors define a particular combination of levels.
     * Combinations are guaranteed to be sorted.
     */
    std::vector<std::vector<Factor> > factors;

    /**
     * Number of cells in each unique combination of factor levels.
     * This has the same length as each inner vector of `factors`.
     * All entries are guaranteed to be positive.
     */
    std::vector<size_t> counts;
};

/**
 * @tparam Factor Factor type.
 * Any type may be used here as long as it is comparable.
 * @tparam Combined Integer type for the combined factor.
 * This should be large enough to hold the number of unique combinations.
 *
 * @param n Number of observations (i.e., cells).
 * @param[in] factors Pointers to arrays of length `n`, each containing a different factor.
 * @param[out] combined Pointer to an array of length `n`, in which the combined factor is to be stored.
 *
 * @return 
 * A `Combinations` object is returned containing the unique combinations of levels observed in `factors`.
 * A combined factor is saved to `combined`, where each entry is an index into the relevant combination of the output `Combinations` object.
 */
template<typename Factor, typename Combined>
Combinations<Factor> combine_factors(size_t n, const std::vector<const Factor*>& factors, Combined* combined) {
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

    std::map<size_t, Combined, decltype(cmp)> mapping(cmp);
    for (size_t i = 0; i < n; ++i) {
        auto mIt = mapping.find(i);
        if (mIt == mapping.end() || !eq(i, mIt->first)) {
            mapping.insert(mIt, std::pair<size_t, Combined>(i, 0));
        }
    }

    // Obtaining the sorted set of unique combinations.
    size_t nfac = factors.size();
    Combinations<Factor> output(nfac);
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
