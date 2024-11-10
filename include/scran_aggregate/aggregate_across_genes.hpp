#ifndef SCRAN_AGGREGATE_AGGREGATE_ACROSS_GENES_HPP
#define SCRAN_AGGREGATE_AGGREGATE_ACROSS_GENES_HPP

#include <algorithm>
#include <vector>
#include <unordered_set>

#include "tatami/tatami.hpp"

/**
 * @file aggregate_across_genes.hpp
 * @brief Aggregate expression values across genes.
 */

namespace scran_aggregate {

/**
 * @brief Options for `aggregate_across_genes()`.
 */
struct AggregateAcrossGenesOptions {
    /**
     * Number of threads to use. 
     * The parallelization scheme is determined by `tatami::parallelize()`.
     */
    int num_threads = 1;

    /**
     * Whether to compute the average expression within each gene set.
     * If the gene set contains weights, a weighted average is computed.
     */
    bool average = false;
};

/**
 * @brief Buffers for `aggregate_across_genes()`.
 * @tparam Sum_ Type of the sum, should be numeric.
 */
template <typename Sum_>
struct AggregateAcrossGenesBuffers {
    /**
     * Vector of length equal to the number of gene sets.
     * Each element is a pointer to an array of length equal to the number of cells,
     * to be filled with the (weighted) sum of expression values for each gene set.
     */
    std::vector<Sum_*> sum;
};

/**
 * @brief Results of `aggregate_across_genes()`.
 * @tparam Sum_ Type of the sum, should be numeric.
 */
template <typename Sum_>
struct AggregateAcrossGenesResults {
    /**
     * Vector of length equal to the number of gene sets.
     * Each inner vector is of length equal to the number of cells.
     * Each entry contains the (weighted) sum of expression values across all genes in the corresponding gene set.
     */
    std::vector<std::vector<Sum_> > sum;
};

/**
 * @cond
 */
namespace aggregate_across_genes_internal {

template<typename Index_, typename Gene_, typename Weight_>
std::vector<Gene_> create_subset(const std::vector<std::tuple<size_t, const Gene_*, const Weight_*> >& gene_sets) {
    std::unordered_set<Gene_> of_interest;
    for (const auto& set : gene_sets) {
        auto set_size = std::get<0>(set);
        auto set_genes = std::get<1>(set);
        of_interest.insert(set_genes, set_genes + set_size);
    }
    std::vector<Index_> subset(of_interest.begin(), of_interest.end());
    std::sort(subset.begin(), subset.end());
    return subset;
}

template<typename Index_>
std::pair<std::vector<Index_>, Index_> create_subset_mapping(const std::vector<Index_>& subset) {
    Index_ offset = 0;
    size_t span = subset.back() - offset + 1;
    std::vector<Index_> mapping(span);
    size_t nsubs = subset.size();
    for (size_t i = 0; i < nsubs; ++i) {
        mapping[subset[i] - offset] = i;
    }
    return std::make_pair(std::move(mapping), offset);
}

template<typename Data_, typename Index_, typename Gene_, typename Weight_, typename Sum_>
void compute_aggregate_by_column(
    const tatami::Matrix<Data_, Index_>& p,
    const std::vector<std::tuple<size_t, const Gene_*, const Weight_*> >& gene_sets,
    const AggregateAcrossGenesBuffers<Sum_>& buffers,
    const AggregateAcrossGenesOptions& options)
{
    // Identifying the subset of rows that actually need to be extracted.
    tatami::VectorPtr<Index_> subset_of_interest = std::make_shared<std::vector<Index_> >(create_subset<Index_>(gene_sets));
    const auto& subset = *subset_of_interest;
    size_t nsubs = subset.size();

    // Creating a mapping back to the gene indices in the subset.
    const size_t num_sets = gene_sets.size();
    std::vector<std::pair<std::vector<Index_>, const Weight_*> > remapping(num_sets);
    if (nsubs) {
        auto sub_mapping = create_subset_mapping(subset);
        const auto& mapping = sub_mapping.first;
        Gene_ offset = sub_mapping.second;

        for (size_t s = 0; s < num_sets; ++s) {
            const auto& set = gene_sets[s];
            auto set_size = std::get<0>(set);
            auto set_genes = std::get<1>(set);

            auto& remapped = remapping[s].first;
            remapped.reserve(set_size);
            for (size_t g = 0; g < set_size; ++g) {
                remapped.push_back(mapping[set_genes[g] - offset]);
            }
            remapping[s].second = std::get<2>(set);
        }
    }

    tatami::parallelize([&](size_t, Index_ start, Index_ length) {
        // We extract as sparse even if it is dense, as it's just
        // easier to index from a dense vector.
        auto ext = tatami::consecutive_extractor<false>(&p, false, start, length, subset_of_interest);
        std::vector<Data_> vbuffer(nsubs);

        for (Index_ x = start, end = start + length; x < end; ++x) {
            auto ptr = ext->fetch(vbuffer.data());
            for (size_t s = 0; s < num_sets; ++s) {
                const auto& set = remapping[s];

                Sum_ value = 0;
                if (set.second) {
                    for (size_t i = 0, send = set.first.size(); i < send; ++i) {
                        value += ptr[set.first[i]] * set.second[i];
                    }
                } else {
                    for (auto ix : set.first) {
                        value += ptr[ix];
                    }
                }

                buffers.sum[s][x] = value;
            }
        }

    }, p.ncol(), options.num_threads);
}

template<typename Data_, typename Index_, typename Gene_, typename Weight_, typename Sum_>
void compute_aggregate_by_row(
    const tatami::Matrix<Data_, Index_>& p,
    const std::vector<std::tuple<size_t, const Gene_*, const Weight_*> >& gene_sets,
    const AggregateAcrossGenesBuffers<Sum_>& buffers,
    const AggregateAcrossGenesOptions& options)
{
    // Identifying the subset of rows that actually need to be extracted.
    auto subset = create_subset<Index_>(gene_sets);
    size_t nsubs = subset.size();
    auto sub_oracle = std::make_shared<tatami::FixedViewOracle<Index_> >(subset.data(), subset.size());

    const size_t num_sets = gene_sets.size();
    std::vector<std::vector<std::pair<size_t, Weight_> > > remapping(nsubs);
    if (nsubs) {
        auto sub_mapping = create_subset_mapping(subset);
        const auto& mapping = sub_mapping.first;
        Gene_ offset = sub_mapping.second;

        for (size_t s = 0; s < num_sets; ++s) {
            const auto& set = gene_sets[s];
            auto set_size = std::get<0>(set);
            auto set_genes = std::get<1>(set);
            auto set_weights = std::get<2>(set);

            if (set_weights) {
                for (size_t g = 0; g < set_size; ++g) {
                    remapping[mapping[set_genes[g] - offset]].emplace_back(s, set_weights[g]); 
                }
            } else {
                for (size_t g = 0; g < set_size; ++g) {
                    remapping[mapping[set_genes[g] - offset]].emplace_back(s, 1);
                }
            }
        }
    }

    Index_ NC = p.ncol();
    if (p.sparse()) {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) {
            auto ext = tatami::new_extractor<true, true>(&p, true, sub_oracle, start, length);
            std::vector<Data_> vbuffer(length);
            std::vector<Index_> ibuffer(length);

            for (const auto& sets : remapping) {
                auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                for (Index_ c = 0; c < range.number; ++c) {
                    auto cell = range.index[c];
                    auto val = range.value[c];
                    for (const auto& s : sets) {
                        buffers.sum[s.first][cell] += val * s.second;
                    }
                }
            }
        }, NC, options.num_threads);

    } else {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) {
            auto ext = tatami::new_extractor<false, true>(&p, true, sub_oracle, start, length);
            std::vector<Data_> vbuffer(length);

            for (const auto& sets : remapping) {
                auto ptr = ext->fetch(vbuffer.data());
                for (Index_ cell = 0; cell < length; ++cell) {
                    auto val = ptr[cell];
                    size_t pos = cell + start;
                    for (const auto& s : sets) {
                        buffers.sum[s.first][pos] += val * s.second;
                    }
                }
            }
        }, NC, options.num_threads);
    }
}

}
/**
 * @endcond
 */

/**
 * Aggregate expression values across gene sets for each cell.
 * This is used to compute a sum/mean of expression values for one or more gene sets/signatures.
 * Each gene in the set can also be weighted, e.g., to account for the strength of regulatory relationships.
 *
 * @tparam Data_ Type of data in the input matrix, should be numeric.
 * @tparam Index_ Integer type of index in the input matrix.
 * @tparam Gene_ Integer type for the indices of genes in each set.
 * @tparam Weight_ Floating-point type for the weights of genes in each set.
 * @tparam Sum_ Floating-point type of the sum.
 *
 * @param input The input matrix where rows are features and columns are cells.
 * @param gene_sets Vector of gene sets.
 * Each tuple corresponds to a set and contains (i) the number of genes in the set,
 * (ii) a pointer to the row indices of the genes in the set, and
 * (iii) a pointer to the weights of the genes in the set.
 * The weight pointer may be NULL, in which case all weights are set to 1.
 * @param[out] buffers Collection of buffers in which to store the aggregate statistics (e.g., sums) for each gene set and cell.
 * @param options Further options.
 */
template<typename Data_, typename Index_, typename Gene_, typename Weight_, typename Sum_>
void aggregate_across_genes(
    const tatami::Matrix<Data_, Index_>& input,
    const std::vector<std::tuple<size_t, const Gene_*, const Weight_*> >& gene_sets,
    const AggregateAcrossGenesBuffers<Sum_>& buffers,
    const AggregateAcrossGenesOptions& options)
{
    if (input.prefer_rows()) {
        aggregate_across_genes_internal::compute_aggregate_by_row(input, gene_sets, buffers, options);
    } else {
        aggregate_across_genes_internal::compute_aggregate_by_column(input, gene_sets, buffers, options);
    }

    if (options.average) {
        size_t nsets = gene_sets.size();
        tatami::parallelize([&](int, size_t start, size_t length) {
            size_t NC = input.ncol();
            for (size_t s = start, end = start + length; s < end; ++s) {
                const auto& set = gene_sets[s];
                auto set_size = std::get<0>(set);

                Sum_ denom = 0;
                auto set_weights = std::get<2>(set);
                if (set_weights) {
                    denom = std::accumulate(set_weights, set_weights + set_size, static_cast<Sum_>(0)); 
                } else {
                    denom = set_size;
                }

                auto current = buffers.sum[s];
                for (size_t c = 0; c < NC; ++c) {
                    current[c] /= denom;
                }
            }
        }, nsets, options.num_threads);
    }
} 

/**
 * Overload of `aggregate_across_genes()` that allocates memory for the results.
 *
 * @tparam Sum_ Floating-point type of the sum.
 * @tparam Data_ Type of data in the input matrix, should be numeric.
 * @tparam Index_ Integer type of index in the input matrix.
 * @tparam Gene_ Integer type for the indices of genes in each set.
 * @tparam Weight_ Floating-point type for the weights of genes in each set.
 *
 * @param input The input matrix where rows are features and columns are cells.
 * @param gene_sets Vector of gene sets.
 * Each tuple corresponds to a set and contains (i) the number of genes in the set,
 * (ii) a pointer to the row indices of the genes in the set, and
 * (iii) a pointer to the weights of the genes in the set.
 * The weight pointer may be NULL, in which case all weights are set to 1.
 * @param options Further options.
 *
 * @return Results of the aggregation.
 */
template<typename Sum_ = double, typename Data_, typename Index_, typename Gene_, typename Weight_>
AggregateAcrossGenesResults<Sum_> aggregate_across_genes(
    const tatami::Matrix<Data_, Index_>& input,
    const std::vector<std::tuple<size_t, const Gene_*, const Weight_*> >& gene_sets,
    const AggregateAcrossGenesOptions& options)
{
    AggregateAcrossGenesResults<Sum_> output;
    AggregateAcrossGenesBuffers<Sum_> buffers;

    size_t NC = input.ncol();
    size_t nsets = gene_sets.size();
    output.sum.resize(nsets);
    buffers.sum.resize(nsets);

    for (size_t s = 0; s < nsets; ++s) {
        output.sum[s].resize(NC);
        buffers.sum[s] = output.sum[s].data();
    }

    aggregate_across_genes(input, gene_sets, buffers, options);
    return output;
} 

}

#endif
