#include "scran_tests/scran_tests.hpp"

#include <map>
#include <random>

#include "scran_aggregate/aggregate_across_cells.hpp"
#include "tatami_stats/tatami_stats.hpp"

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

TEST_P(AggregateAcrossCellsTest, Basics) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    std::vector<int> groupings = create_groupings(dense_row->ncol(), ngroups);

    scran_aggregate::AggregateAcrossCellsOptions opt;
    auto ref = scran_aggregate::aggregate_across_cells(*dense_row, groupings.data(), ngroups, opt);

    auto compare = [&](const auto& other) -> void {
        for (int l = 0; l < ngroups; ++l) {
            scran_tests::compare_almost_equal_containers(ref.sum[l], other.sum[l], {});
            EXPECT_EQ(ref.detected[l], other.detected[l]);
        }
    };

    opt.num_threads = nthreads; 
    if (nthreads != 1) {
        auto res1 = scran_aggregate::aggregate_across_cells(*dense_row, groupings.data(), ngroups, opt);
        compare(res1);
    } else {
        // Doing some cursory checks.
        EXPECT_EQ(ref.sum.size(), ngroups);
        EXPECT_EQ(ref.detected.size(), ngroups);

        std::vector<std::vector<int> > collected_indices(ngroups);
        const int NC = dense_row->ncol();
        for (int c = 0; c < NC; ++c) {
            collected_indices[groupings[c]].push_back(c);
        }

        for (int l = 0; l < ngroups; ++l) {
            auto submat = tatami::make_DelayedSubset(dense_row, std::move(collected_indices[l]), false);
            auto expected_sum = tatami_stats::sum(true, *submat, {});
            scran_tests::compare_almost_equal_containers(expected_sum, ref.sum[l], {});

            tatami::DelayedUnaryIsometricOperation<int, double, int> positive(
                submat, 
                std::make_shared<tatami::DelayedUnaryIsometricGreaterThanScalarHelper<int, double, int, double> >(0.0)
            );
            auto expected_detected = tatami_stats::sum<int>(true, positive, {});
            EXPECT_EQ(expected_detected, ref.detected[l]);
        }
    }

    auto res2 = scran_aggregate::aggregate_across_cells(*sparse_row, groupings.data(), ngroups, opt);
    compare(res2);

    auto res3 = scran_aggregate::aggregate_across_cells(*dense_column, groupings.data(), ngroups, opt);
    compare(res3);

    auto res4 = scran_aggregate::aggregate_across_cells(*sparse_column, groupings.data(), ngroups, opt);
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

/*********************************************/

TEST(AggregateAcrossCells, Skipping) {
    int nr = 88, nc = 126;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulateVectorParameters sparams;
        sparams.density = 0.1;
        sparams.seed = 69;
        return sparams;
    }());
    auto input = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));
    auto grouping = create_groupings(input->ncol(), 2);

    scran_aggregate::AggregateAcrossCellsOptions opt;
    auto ref = scran_aggregate::aggregate_across_cells(*input, grouping.data(), 2, opt);
    EXPECT_EQ(ref.sum.size(), 2);
    EXPECT_EQ(ref.detected.size(), 2);

    // Skipping works correctly when we don't want to compute things.
    opt.compute_sum = false;
    auto partial = scran_aggregate::aggregate_across_cells(*input, grouping.data(), 2, opt);
    EXPECT_EQ(partial.sum.size(), 0);
    EXPECT_EQ(partial.detected.size(), 2);
    
    opt.compute_detected = false;
    auto skipped = scran_aggregate::aggregate_across_cells(*input, grouping.data(), 2, opt);
    EXPECT_EQ(skipped.sum.size(), 0);
    EXPECT_EQ(skipped.detected.size(), 0);
}

/*********************************************/

class AggregateAcrossCellsMedianTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        int nr = 112, nc = 78;
        auto vec = scran_tests::simulate_vector(nr * nc, []{
            scran_tests::SimulateVectorParameters sparams;
            sparams.density = 0.5; // set to 0.5 to encourage a few values above and below for a proper median test.
            sparams.lower = 1;
            sparams.upper = 10;
            return sparams;
        }());

        dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));
        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }
};

TEST_P(AggregateAcrossCellsMedianTest, Basic) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    std::vector<int> groupings = create_groupings(dense_row->ncol(), ngroups);

    scran_aggregate::AggregateAcrossCellsOptions opt;
    opt.compute_median = true;
    auto med = scran_aggregate::aggregate_across_cells(*dense_row, groupings.data(), ngroups, opt);

    auto compare = [&](const auto& other) -> void {
        for (int l = 0; l < ngroups; ++l) {
            EXPECT_EQ(med.sum[l], other.sum[l]);
            EXPECT_EQ(med.detected[l], other.detected[l]);
            EXPECT_EQ(med.median[l], other.median[l]);
        }
    };

    if (nthreads != 1) {
        // Check that the results on parallelization are the same.
        opt.num_threads = nthreads; 
        auto res1 = scran_aggregate::aggregate_across_cells(*dense_row, groupings.data(), ngroups, opt);
        compare(res1);
    } else {
        EXPECT_EQ(med.sum.size(), ngroups);
        EXPECT_EQ(med.detected.size(), ngroups);
        EXPECT_EQ(med.median.size(), ngroups);

        std::vector<std::vector<int> > collected_indices(ngroups);
        const int NC = dense_row->ncol();
        for (int c = 0; c < NC; ++c) {
            collected_indices[groupings[c]].push_back(c);
        }

        for (int l = 0; l < ngroups; ++l) {
            auto submat = tatami::make_DelayedSubset(dense_row, std::move(collected_indices[l]), false);
            auto expected_med = tatami_stats::median(true, *submat, {});
            scran_tests::compare_almost_equal_containers(expected_med, med.median[l], {});

            auto expected_sum = tatami_stats::sum(true, *submat, {});
            scran_tests::compare_almost_equal_containers(expected_sum, med.sum[l], {});

            tatami::DelayedUnaryIsometricOperation<int, double, int> positive(
                submat, 
                std::make_shared<tatami::DelayedUnaryIsometricGreaterThanScalarHelper<int, double, int, double> >(0.0)
            );
            auto expected_detected = tatami_stats::sum<int>(true, positive, {});
            EXPECT_EQ(expected_detected, med.detected[l]);
        }
    }

    auto res2 = scran_aggregate::aggregate_across_cells(*sparse_row, groupings.data(), ngroups, opt);
    compare(res2);

    auto res3 = scran_aggregate::aggregate_across_cells(*dense_column, groupings.data(), ngroups, opt);
    compare(res3);

    auto res4 = scran_aggregate::aggregate_across_cells(*sparse_column, groupings.data(), ngroups, opt);
    compare(res4);
}

INSTANTIATE_TEST_SUITE_P(
    AggregateAcrossCells,
    AggregateAcrossCellsMedianTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of clusters
        ::testing::Values(1, 3) // number of threads
    )
);
