#include "scran_tests/scran_tests.hpp"

#include <map>
#include <random>

#include "scran_aggregate/aggregate_across_genes.hpp"

class AggregateAcrossGenesTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        int nr = 112, nc = 78;
        auto vec = scran_tests::simulate_vector(nr * nc, []{
            scran_tests::SimulateVectorParameters sparams;
            sparams.density = 0.1;
            return sparams;
        }());

        dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));
        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }
};

static std::vector<std::vector<int> > create_gene_sets(int ngenes, int scenario, unsigned long long seed) {
    // Each set is empty.
    if (scenario == 0) {
        return std::vector<std::vector<int> >(19);
    }

    // No sets at all.
    if (scenario == 1) {
        return std::vector<std::vector<int> >();
    }

    // Thread-specific sets.
    if (scenario == 2) {
        const int nsets = 11;
        std::mt19937_64 rng(seed);
        std::vector<std::vector<int> > mock_sets(nsets);

        // Here, we create gene sets where each set only contains a small range of row indices. 
        // This checks the behavior of the parallelized row-major algorithm where each thread processes a separate subset of genes.
        // Some threads will not process any genes for particular sets, in which case they should not allocate temporary memory for those sets.
        // Our aim here is to check the code that skips the memory allocation.
        const int per_set = ngenes / nsets;
        const int remainder = ngenes % nsets;
        for (int s = 0; s < nsets; ++s) {
            const int start = s * per_set + (s < remainder ? s : remainder);
            const int len = per_set + (s < remainder);
            for (int l = 0; l < len; ++l) {
                mock_sets[s].push_back(start + l);
            }
            std::shuffle(mock_sets[s].begin(), mock_sets[s].end(), rng); // shuffling for some variety.
        }
        std::shuffle(mock_sets.begin(), mock_sets.end(), rng); // shuffling for even some variety.
        return mock_sets;
    }

    int nsets, start_gene, end_gene, gene_step;
    double density;
    if (scenario == 3) {
        nsets = 100;
        // Sampling 15% of every tenth gene for each set.
        // This tests that we behave correctly when the subset is not a contiguous block, such that remapping is non-trivial.
        start_gene = 0;
        end_gene = ngenes;
        gene_step = 10;
        density = 0.15;

    } else if (scenario == 4) {
        nsets = 50;
        // Sampling 20% of every third gene for each set.
        // This tests that we behave correctly when the subset is not a contiguous block, such that remapping is non-trivial.
        start_gene = 0;
        end_gene = ngenes;
        gene_step = 3;
        density = 0.2;

    } else if (scenario == 5) { 
        nsets = 79;
        // Ensuring that each gene is represented in at least one set.
        // This checks the special case where all genes are present in the subset, and thus no remapping is required.
        start_gene = 0;
        end_gene = ngenes;
        gene_step = 1;
        density = 0.05;

    } else if (scenario == 6) {
        nsets = 54;
        // Each of the first 25% of genes is represented in at least one set.
        // The aim is to check that we correctly handle contiguous blocks starting at the first gene.
        start_gene = 0;
        end_gene = ngenes / 4;
        gene_step = 1;
        density = 0.1;

    } else {
        nsets = 67;
        // Each of the middle third of genes is represented in at least one set.
        // The aim is to check that we correctly handle contiguous blocks starting after the first gene.
        start_gene = ngenes / 3;
        end_gene = (ngenes * 2) / 3;
        gene_step = 1;
        density = 0.15;
    }

    std::mt19937_64 rng(seed);
    std::vector<std::vector<int> > mock_sets(nsets);

    // Guarantee that each 'gene_step'-th gene in [start_gene, end_gene) is present.
    for (int g = start_gene; g < end_gene; g += gene_step) {
        mock_sets[rng() % nsets].push_back(g);
    }

    // Sprinkling in some more genes for variety.
    // Some care is required to ensure that these remain unique.
    std::uniform_real_distribution runif;
    for (auto& grp : mock_sets) {
        std::unordered_set<int> copy(grp.begin(), grp.end());
        for (int g = start_gene; g < end_gene; g += gene_step) {
            if (runif(rng) < density) {
                copy.insert(g);
            }
        }
        grp.clear();
        grp.insert(grp.end(), copy.begin(), copy.end());
        std::shuffle(grp.begin(), grp.end(), rng);
    }
    return mock_sets;
}

TEST_P(AggregateAcrossGenesTest, Unweighted) {
    auto params = GetParam();
    const auto scenario = std::get<0>(params);
    const auto nthreads = std::get<1>(params);

    const auto mock_sets = create_gene_sets(dense_row->nrow(), scenario, /* seed = */ (scenario + 17) * nthreads);
    const std::size_t nsets = mock_sets.size();

    std::vector<scran_aggregate::AggregateAcrossGenesSet<int, double> > gene_sets;
    gene_sets.reserve(nsets);
    for (const auto& grp : mock_sets) {
        gene_sets.emplace_back(grp.size(), grp.data(), static_cast<double*>(NULL));
    }

    auto compare = [&](const auto& ref, const auto& other) -> void {
        ASSERT_EQ(ref.sum.size(), other.sum.size());
        for (size_t s = 0; s < nsets; ++s) {
            scran_tests::compare_almost_equal_containers(ref.sum[s], other.sum[s], {});
        }
    };

    scran_aggregate::AggregateAcrossGenesOptions opt;
    opt.num_threads = nthreads; 
    auto res1 = scran_aggregate::aggregate_across_genes(*dense_row, gene_sets, opt);
    EXPECT_EQ(res1.sum.size(), nsets);

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
        scran_tests::compare_almost_equal_containers(expected, ave.sum[s], {});
    }
}

TEST_P(AggregateAcrossGenesTest, Weighted) {
    auto params = GetParam();
    const auto scenario = std::get<0>(params);
    const auto nthreads = std::get<1>(params);

    const auto mock_sets = create_gene_sets(dense_row->nrow(), scenario, /* seed = */ (scenario + 13) * nthreads);
    const std::size_t nsets = mock_sets.size();

    std::vector<std::vector<double> > weights(nsets);
    {
        std::mt19937_64 rng((scenario + 17) * nthreads);
        std::uniform_real_distribution runif;
        for (size_t s = 0; s < nsets; ++s) {
            auto& wt = weights[s];
            for ([[maybe_unused]] auto g : mock_sets[s]) {
                wt.push_back(runif(rng));
            }
        }
    }

    std::vector<scran_aggregate::AggregateAcrossGenesSet<int, double> > gene_sets;
    gene_sets.reserve(nsets);
    for (size_t s = 0; s < nsets; ++s) {
        const auto& grp = mock_sets[s];
        gene_sets.emplace_back(grp.size(), grp.data(), weights[s].data());
    }

    auto compare = [&](const auto& ref, const auto& other) -> void {
        ASSERT_EQ(ref.sum.size(), other.sum.size());
        for (size_t s = 0; s < nsets; ++s) {
            scran_tests::compare_almost_equal_containers(ref.sum[s], other.sum[s], {});
        }
    };

    scran_aggregate::AggregateAcrossGenesOptions opt;
    opt.num_threads = nthreads; 
    auto res1 = scran_aggregate::aggregate_across_genes(*dense_row, gene_sets, opt);
    EXPECT_EQ(res1.sum.size(), nsets);

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
        scran_tests::compare_almost_equal_containers(expected, ave.sum[s], {});
    }
}

INSTANTIATE_TEST_SUITE_P(
    AggregateAcrossGenes,
    AggregateAcrossGenesTest,
    ::testing::Combine(
        ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7), // scenarios
        ::testing::Values(1, 3) // number of threads
    )
);

TEST(AggregateAcrossGenes, OutOfRange) {
    int nr = 11, nc = 78;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulateVectorParameters sparams;
        sparams.density = 0.1;
        return sparams;
    }());

    tatami::DenseRowMatrix<double, int> mat(nr, nc, std::move(vec));
    std::vector<int> example { 1, 10, 100 };
    std::vector<scran_aggregate::AggregateAcrossGenesSet<int, double> > gene_sets;
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
