#include "scran_tests/scran_tests.hpp"

#include <map>
#include <random>

#include "scran_aggregate/aggregate_across_cells.hpp"

static std::vector<int> create_groupings(size_t n, int ngroups) {
    std::vector<int> groupings(n);
    for (size_t g = 0; g < groupings.size(); ++g) {
        groupings[g] = g % ngroups;
    }
    return groupings;
}

/*********************************************/

class AggregateAcrossCellsTest : public ::testing::TestWithParam<std::tuple<int, int> > {
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

TEST_P(AggregateAcrossCellsTest, Basics) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    std::vector<int> groupings = create_groupings(dense_row->ncol(), ngroups);

    scran_aggregate::AggregateAcrossCellsOptions opt;
    auto ref = scran_aggregate::aggregate_across_cells(*dense_row, groupings.data(), opt);

    auto compare = [&](const auto& other) -> void {
        for (int l = 0; l < ngroups; ++l) {
            EXPECT_EQ(ref.sums[l], other.sums[l]);
            EXPECT_EQ(ref.detected[l], other.detected[l]);
        }
    };

    opt.num_threads = nthreads; 
    if (nthreads != 1) {
        auto res1 = scran_aggregate::aggregate_across_cells(*dense_row, groupings.data(), opt);
        compare(res1);
    } else {
        // Doing some cursory checks.
        EXPECT_EQ(ref.sums.size(), ngroups);
        EXPECT_EQ(ref.detected.size(), ngroups);

        auto NR = dense_row->nrow();
        for (int l = 0; l < ngroups; ++l) {
            EXPECT_EQ(ref.sums[l].size(), NR);
            EXPECT_EQ(ref.detected[l].size(), NR);
            EXPECT_NE(std::accumulate(ref.sums[l].begin(), ref.sums[l].end(), 0.0), 0); // check that something is present, at least.
            EXPECT_NE(std::accumulate(ref.detected[l].begin(), ref.detected[l].end(), 0), 0);
        }
    }

    auto res2 = scran_aggregate::aggregate_across_cells(*sparse_row, groupings.data(), opt);
    compare(res2);

    auto res3 = scran_aggregate::aggregate_across_cells(*dense_column, groupings.data(), opt);
    compare(res3);

    auto res4 = scran_aggregate::aggregate_across_cells(*sparse_column, groupings.data(), opt);
    compare(res4);
}

INSTANTIATE_TEST_SUITE_P(
    AggregateAcrossCells,
    AggregateAcrossCellsTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of clusters
        ::testing::Values(1, 3) // number of threads
    )
);

TEST(AggregateAcrossCells, Skipping) {
    int nr = 88, nc = 126;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulationParameters sparams;
        sparams.density = 0.1;
        sparams.seed = 69;
        return sparams;
    }());
    auto input = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));
    auto grouping = create_groupings(input->ncol(), 2);

    scran_aggregate::AggregateAcrossCellsOptions opt;
    auto ref = scran_aggregate::aggregate_across_cells(*input, grouping.data(), opt);
    EXPECT_EQ(ref.sums.size(), 2);
    EXPECT_EQ(ref.detected.size(), 2);

    // Skipping works correctly when we don't want to compute things.
    opt.compute_sums = false;
    auto partial = scran_aggregate::aggregate_across_cells(*input, grouping.data(), opt);
    EXPECT_EQ(partial.sums.size(), 0);
    EXPECT_EQ(partial.detected.size(), 2);
    
    opt.compute_detected = false;
    auto skipped = scran_aggregate::aggregate_across_cells(*input, grouping.data(), opt);
    EXPECT_EQ(skipped.sums.size(), 0);
    EXPECT_EQ(skipped.detected.size(), 0);
}
