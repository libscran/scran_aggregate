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
 * Any type may be used here as long as it implements the comparison operators.
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

/**
 * This function is a variation of `combine_factors()` that considers unobserved combinations of factor levels.
 *
 * @tparam Factor_ Factor type.
 * Any type may be used here as long as it is comparable.
 * @tparam Number_ Integer type for the number of levels in each factor.
 * @tparam Combined_ Integer type for the combined factor.
 * This should be large enough to hold the number of unique (possibly unused) combinations.
 *
 * @param n Number of observations (i.e., cells).
 * @param[in] factors Vector of pairs, each of which corresponds to a factor.
 * The first element of the pair is a pointer to an array of length `n`, containing the factor level for each observation.
 * The second element is the total number of levels for this factor, which may be greater than the largeset observed level.
 * @param[out] combined Pointer to an array of length `n` in which the combined factor is to be stored.
 * On output, each entry determines the corresponding observation's combination of levels by indexing into the inner vectors of the returned object;
 * see the argument of the same name in `combine_factors()` for more details.
 *
 * @return 
 * Vector of vectors containing each unique combinations of factor levels.
 * This has the same structure as the output of `combine_factors()`,
 * with the only difference being that unobserved combinations are also reported.
 */
template<typename Factor_, typename Number_, typename Combined_>
std::vector<std::vector<Factor_> > combine_factors_unused(size_t n, const std::vector<std::pair<const Factor_*, Number_> >& factors, Combined_* combined) {
    size_t nfac = factors.size();
    std::vector<std::vector<Factor_> > output(nfac);

    if (nfac > 1) {
        // We iterate from back to front, where the first factor is the slowest changing.
        std::copy_n(factors[nfac - 1].first, n, combined); 
        Combined_ mult = factors[nfac - 1].second;
        for (size_t f = nfac - 1; f > 0; --f) {
            const auto& finfo = factors[f - 1];
            auto ff = finfo.first;
            for (size_t i = 0; i < n; ++i) {
                combined[i] += mult * ff[i];
            }
            mult *= finfo.second;
        }

        auto ncombos = mult;
        Combined_ outer_repeats = mult;
        Combined_ inner_repeats = 1;
        for (size_t f = nfac; f > 0; --f) {
            auto& out = output[f - 1];
            out.reserve(ncombos);

            const auto& finfo = factors[f - 1];
            size_t initial_size = inner_repeats * finfo.second;
            out.resize(initial_size);

            if (inner_repeats == 1) {
                std::iota(out.begin(), out.end(), static_cast<Combined_>(0));
            } else {
                auto oIt = out.begin();
                for (Number_ l = 0; l < finfo.second; ++l) {
                    std::fill_n(oIt, inner_repeats, l);
                    oIt += inner_repeats;
                }
            }
            inner_repeats = initial_size;

            outer_repeats /= finfo.second;
            for (Combined_ r = 1; r < outer_repeats; ++r) {
                out.insert(out.end(), out.begin(), out.begin() + initial_size);
            }
        }

    } else if (nfac == 1) {
        output[0].resize(factors[0].second);
        std::iota(output[0].begin(), output[0].end(), static_cast<Combined_>(0));
        std::copy_n(factors[0].first, n, combined);

    } else {
        std::fill_n(combined, n, 0);
    }

    return output;
}

}

#endif
