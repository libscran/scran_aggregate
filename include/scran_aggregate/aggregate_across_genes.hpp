#ifndef SCRAN_AGGREGATE_AGGREGATE_ACROSS_GENES_HPP
#define SCRAN_AGGREGATE_AGGREGATE_ACROSS_GENES_HPP

#include <algorithm>
#include <vector>
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
    std::vector<Detected_*> sum;
};

/**
 * @brief Results of `aggregate_across_genes()`.
 * @tparam Sum_ Type of the sum, should be numeric.
 */
template <typename Sum_, typename Detected_>
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
namespace internal {

template<typename Gene_, typename Weight_>
std::unordered_set<Gene_> create_gene_universe(const std::vector<std::tuple<size_t, const Gene_*, const Weight_*> >& gene_sets) {
    std::unordered_set<Gene_> of_interest;
    for (const auto& set : gene_sets) {
        auto set_size = std::get<0>(set);
        auto set_genes = std::get<1>(set);
        of_interest.insert(set_genes, set_genes + set_size);
    }
    return of_interest;
}

template<typename Index_>
std::pair<std::vector<Index_>, Index_> create_subset_mapping(const std::vector<Index_>& subset) {
    Index_ offset = 0;
    size_t span = subset.back() - offset + 1;
    std::vector<Index_> mapping(span);
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
    auto of_interest = create_gene_universe(gene_sets);
    tatami::VectorPtr<Index_> subset_of_interest = std::make_shared<std::vector<Index_> >(of_interest.begin(), of_interest.end());
    const auto& subset = *subset_of_interest;
    size_t nsubs = subset.size();

    // Creating a mapping back to the gene indices in the subset.
    const size_t num_sets = gene_sets.size();
    std::vector<std::pair<std::vector<Index_>, const Weight_*> > remapping(num_sets);
    if (nsubs) {
        auto sub_mapping = create_mapping(subset);
        const auto& mapping = sub_mapping.first;
        Gene_ offset = mapping.second;

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
        auto ext = tatami::consecutive_extractor<true>(&p, false, start, length, subset_of_interest);
        auto NR = p.nrow();
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
    auto of_interest = create_gene_universe(gene_sets);
    std::vector<Index_> subset(of_interest.begin(), of_interest.end());
    size_t nsubs = subset.size();
    auto sub_oracle = std::make_shared<tatami::FixedViewOracle<Index_> >(subset.data(), subset.size());

    std::vector<std::vector<std::pair<size_t, Weight_> > > remapping(nsubs);
    if (nsubs) {
        auto sub_mapping = create_mapping(subset);
        const auto& mapping = sub_mapping.first;
        Gene_ offset = mapping.second;

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
            for (Index_ x = start, end = start + length; x < end; ++x) {
                auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                const auto& sets = remapping[x]; 
                for (Index_ c = 0; c < range.number; ++c) {
                    auto cell = range.index[c];
                    auto val = range.value[c];
                    for (const auto& s : sets) {
                        buffers.sums[s.first][cell] += val * s.second;
                    }
                }
            }
        }, NC, options.num_threads);

    } else {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) {
            auto ext = tatami::new_extractor<false, true>(&p, true, sub_oracle, start, length);
            std::vector<Data_> vbuffer(length);
            for (Index_ x = start, end = start + length; x < end; ++x) {
                auto ptr = ext->fetch(vbuffer.data());
                const auto& sets = remapping[x]; 
                for (Index_ cell = 0; cell < NC; ++cell) {
                    auto val = ptr[cell];
                    for (const auto& s : sets) {
                        buffers.sums[s.first][cell] += val * s.second;
                    }
                }
            }
        }, NC, options.num_threads);
    }

}


/**
 * @endcond
 */

}

}

#endif
