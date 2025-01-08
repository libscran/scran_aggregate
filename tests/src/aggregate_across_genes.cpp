#include <gtest/gtest.h>

#include "scran_aggregate/aggregate_across_genes.hpp"
#include "scran_tests/scran_tests.hpp"
#include <map>
#include <random>

static std::vector<std::vector<int> > create_sets(size_t nsets, int ngenes, size_t seed) {
    std::vector<std::vector<int> > groupings(nsets);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution runif;
    for (auto& grp : groupings) {
        for (int g = 0; g < ngenes; ++g) {
            if (runif(rng) < 0.15) {
                grp.push_back(g);
            }
        }
    }
    return groupings;
}

/*********************************************/

class AggregateAcrossGenesTest : public ::testing::TestWithParam<int> {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        int nr = 112, nc = 78;
        auto vec = scran_tests::simulate_vector(nr * nc, []{
            scran_tests::SimulationParameters sparams;
            sparams.density = 0.1;
            return sparams;
        }());

        dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));
        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }
};

TEST_P(AggregateAcrossGenesTest, Unweighted) {
    auto nthreads = GetParam();

    const size_t nsets = 100;
    int ngenes = dense_row->nrow();
    auto mock_sets = create_sets(nsets, ngenes, /* seed = */ nsets * nthreads);

    std::vector<std::tuple<size_t, const int*, const double*> > gene_sets;
    gene_sets.reserve(mock_sets.size());
    for (const auto& grp : mock_sets) {
        gene_sets.emplace_back(grp.size(), grp.data(), static_cast<double*>(NULL));
    }

    auto compare = [&](const auto& ref, const auto& other) -> void {
        for (size_t s = 0; s < nsets; ++s) {
            EXPECT_EQ(ref.sum[s], other.sum[s]);
        }
    };

    scran_aggregate::AggregateAcrossGenesOptions opt;
    opt.num_threads = nthreads; 
    auto res1 = scran_aggregate::aggregate_across_genes(*dense_row, gene_sets, opt);

    if (nthreads > 1) {
        auto copy = opt;
        copy.num_threads = 1;
        auto ref = scran_aggregate::aggregate_across_genes(*dense_row, gene_sets, copy);
        compare(res1, ref);
    }

    auto res2 = scran_aggregate::aggregate_across_genes(*sparse_row, gene_sets, opt);
    compare(res1, res2);

    auto res3 = scran_aggregate::aggregate_across_genes(*dense_column, gene_sets, opt);
    compare(res1, res3);

    auto res4 = scran_aggregate::aggregate_across_genes(*sparse_column, gene_sets, opt);
    compare(res1, res4);

    // Checking that the average works.
    opt.average = true;
    auto ave = scran_aggregate::aggregate_across_genes(*sparse_column, gene_sets, opt);
    for (size_t s = 0; s < nsets; ++s) {
        auto expected = res1.sum[s];
        for (auto& x : expected) { x /= mock_sets[s].size(); }
        EXPECT_EQ(expected, ave.sum[s]);
    }
}

TEST_P(AggregateAcrossGenesTest, Weighted) {
    auto nthreads = GetParam();

    size_t nsets = 50;
    int ngenes = dense_row->nrow();
    std::vector<std::vector<int> > mock_sets(nsets);
    std::vector<std::vector<double> > weights(nsets);
    {
        std::mt19937_64 rng(nsets * nthreads + 17);
        std::uniform_real_distribution runif;
        for (size_t s = 0; s < nsets; ++s) {
            auto& grp = mock_sets[s];
            auto& wt = weights[s];
            for (int g = 0; g < ngenes; ++g) {
                if (runif(rng) < 0.15) {
                    grp.push_back(g);
                    wt.push_back(runif(rng));
                }
            }
        }
    }

    std::vector<std::tuple<size_t, const int*, const double*> > gene_sets;
    gene_sets.reserve(nsets);
    for (size_t s = 0; s < nsets; ++s) {
        const auto& grp = mock_sets[s];
        gene_sets.emplace_back(grp.size(), grp.data(), weights[s].data());
    }

    auto compare = [&](const auto& ref, const auto& other) -> void {
        for (size_t s = 0; s < nsets; ++s) {
            EXPECT_EQ(ref.sum[s], other.sum[s]);
        }
    };

    scran_aggregate::AggregateAcrossGenesOptions opt;
    opt.num_threads = nthreads; 
    auto res1 = scran_aggregate::aggregate_across_genes(*dense_row, gene_sets, opt);

    if (nthreads > 1) {
        auto copy = opt;
        copy.num_threads = 1;
        auto ref = scran_aggregate::aggregate_across_genes(*dense_row, gene_sets, copy);
        compare(res1, ref);
    }

    auto res2 = scran_aggregate::aggregate_across_genes(*sparse_row, gene_sets, opt);
    compare(res1, res2);

    auto res3 = scran_aggregate::aggregate_across_genes(*dense_column, gene_sets, opt);
    compare(res1, res3);

    auto res4 = scran_aggregate::aggregate_across_genes(*sparse_column, gene_sets, opt);
    compare(res1, res4);

    // Checking that the average works.
    opt.average = true;
    auto ave = scran_aggregate::aggregate_across_genes(*sparse_column, gene_sets, opt);
    for (size_t s = 0; s < nsets; ++s) {
        auto expected = res1.sum[s];
        double denom = std::accumulate(weights[s].begin(), weights[s].end(), 0.0);
        for (auto& x : expected) { x /= denom; }
        EXPECT_EQ(expected, ave.sum[s]);
    }
}

TEST_P(AggregateAcrossGenesTest, Empty) {
    auto nthreads = GetParam();
    std::vector<std::tuple<size_t, const int*, const double*> > gene_sets;

    scran_aggregate::AggregateAcrossGenesOptions opt;
    opt.num_threads = nthreads; 
    auto res1 = scran_aggregate::aggregate_across_genes(*dense_row, gene_sets, opt);
    EXPECT_EQ(res1.sum.size(), 0);

    auto res2 = scran_aggregate::aggregate_across_genes(*sparse_row, gene_sets, opt);
    EXPECT_EQ(res2.sum.size(), 0);

    auto res3 = scran_aggregate::aggregate_across_genes(*dense_column, gene_sets, opt);
    EXPECT_EQ(res3.sum.size(), 0);

    auto res4 = scran_aggregate::aggregate_across_genes(*sparse_column, gene_sets, opt);
    EXPECT_EQ(res4.sum.size(), 0);
}

INSTANTIATE_TEST_SUITE_P(
    AggregateAcrossGenes,
    AggregateAcrossGenesTest,
    ::testing::Values(1, 3) // number of threads
);

TEST(AggregateAcrossGenes, OutOfRange) {
    int nr = 11, nc = 78;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulationParameters sparams;
        sparams.density = 0.1;
        return sparams;
    }());

    tatami::DenseRowMatrix<double, int> mat(nr, nc, std::move(vec));
    std::vector<int> example { 1, 10, 100 };
    std::vector<std::tuple<size_t, const int*, const double*> > gene_sets;
    gene_sets.emplace_back(3, example.data(), static_cast<double*>(NULL));

    scran_aggregate::AggregateAcrossGenesOptions opt;
    scran_tests::expect_error([&]() {
        scran_aggregate::aggregate_across_genes(mat, gene_sets, opt);
    }, "out of range");

    // Also fails if there are negative values.
    example[0] = -1;
    scran_tests::expect_error([&]() {
        scran_aggregate::aggregate_across_genes(mat, gene_sets, opt);
    }, "out of range");
}

TEST(AggregateAcrossGenes, DirtyBuffers) {
    int nr = 110, nc = 78;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulationParameters sparams;
        sparams.density = 0.1;
        return sparams;
    }());

    tatami::DenseRowMatrix<double, int> dense_row(nr, nc, std::move(vec));
    auto sparse_column = tatami::convert_to_compressed_sparse(&dense_row, false);

    // Setting up the gene sets.
    size_t nsets = 100;
    auto mock_sets = create_sets(nsets, nr, /* seed = */ 99);

    std::vector<std::tuple<size_t, const int*, const double*> > gene_sets;
    gene_sets.reserve(mock_sets.size());
    for (const auto& grp : mock_sets) {
        gene_sets.emplace_back(grp.size(), grp.data(), static_cast<double*>(NULL));
    }

    // Setting up some buffers.
    scran_aggregate::AggregateAcrossGenesResults<double> store;
    scran_aggregate::AggregateAcrossGenesBuffers<double> buffers;
    store.sum.resize(nsets);
    buffers.sum.resize(nsets);
    for (size_t s = 0; s < nsets; ++s) {
        store.sum[s].resize(nc, -1); // using -1 to simulate uninitialized values.
        buffers.sum[s] = store.sum[s].data();
    }

    // Checking that the dirtiness doesn't affect the results.
    scran_aggregate::AggregateAcrossGenesOptions opt;
    auto ref = scran_aggregate::aggregate_across_genes(dense_row, gene_sets, opt);
    scran_aggregate::aggregate_across_genes(dense_row, gene_sets, buffers, opt);
    for (size_t s = 0; s < nsets; ++s) {
        EXPECT_EQ(ref.sum[s], store.sum[s]);
    }

    for (size_t s = 0; s < nsets; ++s) {
        std::fill_n(store.sum[s].data(), nc, -1); // dirtying it again for another run.
    }
    scran_aggregate::aggregate_across_genes(*sparse_column, gene_sets, buffers, opt);
    for (size_t s = 0; s < nsets; ++s) {
        EXPECT_EQ(ref.sum[s], store.sum[s]);
    }
}
