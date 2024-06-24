#ifndef SCRAN_AGGREGATE_ACROSS_CELLS_COMPUTE_HPP
#define SCRAN_AGGREGATE_ACROSS_CELLS_COMPUTE_HPP

#include <algorithm>
#include <vector>
#include "tatami/tatami.hpp"

/**
 * @file compute.hpp
 *
 * @brief Aggregate expression values across cells.
 */

namespace aggregate_across_cells {

struct Options {
    /**
     * Whether to compute the sum within each factor level.
     * This option only affects the `compute()` overload where a `Results` object is returned.
     */
    bool compute_sums = true;

    /**
     * Whether to compute the number of detected cells within each factor level.
     * This option only affects the `compute()` overload where a `Results` object is returned.
     */
    bool compute_detected = true;

    /**
     * Number of threads to use. 
     * See `set_num_threads()`.
     */
    int num_threads = 1;
};

/**
 * @cond
 */
namespace internal {

template<bool sparse_, typename Data_, typename Index_, typename Factor_, typename Sum_, typename Detected_>
void compute_by_row(const tatami::Matrix<Data_, Index_>* p, const Factor_* factor, std::vector<Sum_*>& sums, std::vector<Detected_*>& detected, const Options& options) {
    tatami::Options opt;
    opt.sparse_ordered_index = false;

    tatami::parallelize([&](size_t, Index_ s, Index_ l) {
        auto ext = tatami::consecutive_extractor<sparse_>(p, true, s, l, opt);
        std::vector<Sum_> tmp_sums(sums.size());
        std::vector<Detected_> tmp_detected(detected.size());

        auto NC = p->ncol();
        std::vector<Data_> vbuffer(NC);
        typename std::conditional<sparse_, std::vector<Index_>, Index_>::type ibuffer(NC);

        for (Index_ x = s, end = s + l; x < end; ++x) {
            auto row = [&]() {
                if constexpr(sparse_) {
                    return ext->fetch(vbuffer.data(), ibuffer.data());
                } else {
                    return ext->fetch(vbuffer.data());
                }
            }();

            if (sums.size()) {
                std::fill(tmp_sums.begin(), tmp_sums.end(), 0);

                if constexpr(sparse_) {
                    for (Index_ j = 0; j < row.number; ++j) {
                        tmp_sums[factor[row.index[j]]] += row.value[j];
                    }
                } else {
                    for (Index_ j = 0; j < NC; ++j) {
                        tmp_sums[factor[j]] += row[j];
                    }
                }

                // Computing before transferring for more cache-friendliness.
                for (Index_ l = 0; l < tmp_sums.size(); ++l) {
                    sums[l][x] = tmp_sums[l];
                }
            }

            if (detected.size()) {
                std::fill(tmp_detected.begin(), tmp_detected.end(), 0);

                if constexpr(sparse_) {
                    for (Index_ j = 0; j < row.number; ++j) {
                        tmp_detected[factor[row.index[j]]] += (row.value[j] > 0);
                    }
                } else {
                    for (Index_ j = 0; j < NC; ++j) {
                        tmp_detected[factor[j]] += (row[j] > 0);
                    }
                }

                for (Index_ l = 0; l < tmp_detected.size(); ++l) {
                    detected[l][x] = tmp_detected[l];
                }
            }
        }
    }, p->nrow(), options.num_threads);
}

template<bool sparse_, typename Data_, typename Index_, typename Factor_, typename Sum_, typename Detected_>
void compute_by_column(const tatami::Matrix<Data_, Index_>* p, const Factor_* factor, std::vector<Sum_*>& sums, std::vector<Detected_*>& detected, const Options& options) {
    tatami::Options opt;
    opt.sparse_ordered_index = false;

    tatami::parallelize([&](size_t, Index_ s, Index_ l) {
        auto NC = p->ncol();
        auto ext = tatami::consecutive_extractor<sparse_>(p, false, 0, NC, s, l, opt);
        std::vector<Data_> vbuffer(l);
        typename std::conditional<sparse_, std::vector<Index_>, Index_>::type ibuffer(l);

        for (Index_ x = 0; x < NC; ++x) {
            auto current = factor[x];

            if constexpr(sparse_) {
                auto col = ext->fetch(vbuffer.data(), ibuffer.data());
                if (sums.size()) {
                    auto& cursum = sums[current];
                    for (Index_ i = 0; i < col.number; ++i) {
                        cursum[col.index[i]] += col.value[i];
                    }
                }

                if (detected.size()) {
                    auto& curdetected = detected[current];
                    for (Index_ i = 0; i < col.number; ++i) {
                        curdetected[col.index[i]] += (col.value[i] > 0);
                    }
                }

            } else {
                auto col = ext->fetch(vbuffer.data());

                if (sums.size()) {
                    auto cursum = sums[current] + s;
                    for (Index_ i = 0; i < l; ++i) {
                        cursum[i] += col[i];
                    }
                }

                if (detected.size()) {
                    auto curdetected = detected[current] + s;
                    for (Index_ i = 0; i < l; ++i) {
                        curdetected[i] += (col[i] > 0);
                    }
                }
            }
        }
    }, p->nrow(), options.num_threads);
}

}
/**
 * @endcond
 */

/**
 * Compute the per-gene sum and number of cells with detected expression for each level of a grouping factor.
 * This is typically done for the creation of pseudo-bulk expression profiles for cluster/sample combinations.
 * Expression values are generally expected to be counts, though the same code can be trivially re-used to compute the average log-expression.
 * We can also report the number of cells with detected (i.e., positive) expression values in each grouping.
 *
 * @tparam Data_ Type of data in the input matrix, should be numeric.
 * @tparam Index_ Integer type of index in the input matrix.
 * @tparam Factor_ Integer type of the factor.
 * @tparam Sum_ Type of the sum, usually the same as `Data`.
 * @tparam Detected_ Type for the number of detected cells, usually integer.
 *
 * @param input The input matrix where rows are features and columns are cells.
 * @param[in] factor Pointer to an array of length equal to the number of columns of `input`,
 * containing the factor level for each cell.
 * All levels should be integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique levels.
 * @param[out] sums Vector of length \f$N\f$ (see `factor`),
 * containing pointers to arrays of length equal to the number of columns of `input`.
 * These will be filled with the summed expression across all cells in the corresponding level for each gene.
 * Alternatively, if the vector is of length 0, no sums will be computed.
 * @param[out] detected Vector of length \f$N\f$ (see `factor`),
 * containing pointers to arrays of length equal to the number of columns of `input`.
 * These will be filled with the number of cells with detected expression in the corresponding level for each gene.
 * Alternatively, if the vector is of length 0, no numbers will be computed.
 * @param options Further options.
 */
template<typename Data_, typename Index_, typename Factor_, typename Sum_, typename Detected_>
void compute(const tatami::Matrix<Data_, Index_>* input, const Factor_* factor, std::vector<Sum_*> sums, std::vector<Detected_*> detected, const Options& options) {
    if (input->prefer_rows()) {
        if (input->sparse()) {
            internal::compute_by_row<true>(input, factor, sums, detected, options);
        } else {
            internal::compute_by_row<false>(input, factor, sums, detected, options);
        }
    } else {
        if (input->sparse()) {
            internal::compute_by_column<true>(input, factor, sums, detected, options);
        } else {
            internal::compute_by_column<false>(input, factor, sums, detected, options);
        }
    }
} 

/**
 * @brief Container for the aggregation results.
 * @tparam Sum_ Type of the sum, should be numeric.
 * @tparam Detected_ Type for the number of detected cells, usually integer.
 */
template <typename Sum_, typename Detected_>
struct Results {
    /**
     * Vector of length equal to the number of factor levels.
     * Each inner vector is of length equal to the number of genes.
     * Each entry contains the summed expression across all cells in the corresponding level for the corresponding gene.
     *
     * If `Options::compute_sums = false`, this vector is empty.
     */
    std::vector<std::vector<Sum_> > sums;

    /**
     * Vector of length equal to the number of factor levels.
     * Each inner vector is of length equal to the number of genes.
     * Each entry contains the number of cells in the corresponding level with detected expression for the corresponding gene.
     *
     * If `Options::compute_detected = false`, this vector is empty.
     */
    std::vector<std::vector<Detected_> > detected;
};

/**
 * @tparam Sum_ Type of the sum, should be numeric.
 * @tparam Detected_ Type for the number of detected cells, usually integer.
 * @tparam Data_ Type of data in the input matrix, should be numeric.
 * @tparam Index_ Integer type of index in the input matrix.
 * @tparam Factor_ Integer type of the factor.
 *
 * @param input The input matrix where rows are features and columns are cells.
 * @param[in] factor Pointer to an array of length equal to the number of columns of `input`,
 * containing the factor level for each cell.
 * All levels should be integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique levels.
 * @param options Further options.
 *
 * @return A `Results` object is returned, where the available statistics depend on `Options`.
 */
template<typename Sum_ = double, typename Detected_ = int, typename Data_, typename Index_, typename Factor_>
Results<Sum_, Detected_> compute(const tatami::Matrix<Data_, Index_>* input, const Factor_* factor, const Options& options) {
    size_t NC = input->ncol();
    size_t nlevels = (NC ? *std::max_element(factor, factor + NC) + 1 : 0);
    size_t ngenes = input->nrow();

    Results<Sum_, Detected_> output;
    std::vector<Sum_*> sumptr;
    std::vector<Detected_*> detptr;

    if (options.compute_sums) {
        output.sums.resize(nlevels, std::vector<Sum_>(ngenes));
        sumptr.resize(nlevels);
        for (size_t l = 0; l < nlevels; ++l) {
            sumptr[l] = output.sums[l].data();
        }
    }

    if (options.compute_detected) {
        output.detected.resize(nlevels, std::vector<Detected_>(ngenes));
        detptr.resize(nlevels);
        for (size_t l = 0; l < nlevels; ++l) {
            detptr[l] = output.detected[l].data();
        }
    }

    compute(input, factor, std::move(sumptr), std::move(detptr), options);
    return output;
} 

}

#endif
