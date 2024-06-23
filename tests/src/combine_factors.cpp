#include <gtest/gtest.h>

#include <random>

#include "aggregate_across_cells/combine_factors.hpp"

template<typename Factor>
std::pair<aggregate_across_cells::Combinations<Factor>, std::vector<int> > test_combine_factors(size_t n, const std::vector<const Factor*>& factors) {
    std::vector<int> combined(n);
    auto levels = aggregate_across_cells::combine_factors(n, factors, combined.data());
    return std::make_pair(std::move(levels), std::move(combined));
}

TEST(CombineFactors, Simple) {
    // Simple sorted case.
    {
        std::vector<int> stuff{0, 0, 1, 1, 1, 2, 2, 2, 2 };
        auto combined = test_combine_factors(stuff.size(), std::vector<const int*>{stuff.data()});
        EXPECT_EQ(combined.second, stuff);

        EXPECT_EQ(combined.first.factors.size(), 1);
        std::vector<int> expected { 0, 1, 2 };
        EXPECT_EQ(combined.first.factors[0], expected);

        std::vector<size_t> counts { 2, 3, 4 };
        EXPECT_EQ(combined.first.counts, counts);
    }

    // Testing the unsorted case.
    {
        std::vector<int> stuff{ 1, 0, 1, 2, 1, 0, 2, 3, 2 };
        auto combined = test_combine_factors(stuff.size(), std::vector<const int*>{stuff.data()});
        EXPECT_EQ(combined.second, stuff);

        EXPECT_EQ(combined.first.factors.size(), 1);
        std::vector<int> expected { 0, 1, 2, 3 };
        EXPECT_EQ(combined.first.factors[0], expected);

        std::vector<size_t> counts { 2, 3, 3, 1 };
        EXPECT_EQ(combined.first.counts, counts);
    }

    // Non-consecutive still works.
    {
        std::vector<int> stuff{ 1, 3, 5, 7, 9 };
        auto combined = test_combine_factors(stuff.size(), std::vector<const int*>{stuff.data()});
        std::vector<int> expected { 0, 1, 2, 3, 4 };
        EXPECT_EQ(combined.second, expected);

        EXPECT_EQ(combined.first.factors.size(), 1);
        std::vector<int> levels { 1, 3, 5, 7, 9 };
        EXPECT_EQ(combined.first.factors[0], levels);

        std::vector<size_t> counts { 1, 1, 1, 1, 1 };
        EXPECT_EQ(combined.first.counts, counts);
    }
}

TEST(CombineFactors, Multiple) {
    {
        std::vector<int> stuff1{ 0, 0, 1, 1, 1, 2, 2, 2, 2 };
        std::vector<int> stuff2{ 0, 1, 2, 0, 1, 2, 0, 1, 2 };
        auto combined = test_combine_factors(stuff1.size(), std::vector<const int*>{stuff1.data(), stuff2.data()});

        std::vector<int> expected { 0, 1, 4, 2, 3, 7, 5, 6, 7 };
        EXPECT_EQ(combined.second, expected);

        EXPECT_EQ(combined.first.factors.size(), 2);
        std::vector<int> levels1 { 0, 0, 1, 1, 1, 2, 2, 2 };
        std::vector<int> levels2 { 0, 1, 0, 1, 2, 0, 1, 2 };
        EXPECT_EQ(combined.first.factors[0], levels1);
        EXPECT_EQ(combined.first.factors[1], levels2);

        std::vector<size_t> counts { 1, 1, 1, 1, 1, 1, 1, 2 };
        EXPECT_EQ(combined.first.counts, counts);
    }

    {
        std::map<std::pair<int, int>, std::vector<int> > collected;
        std::vector<int> stuff1;
        std::vector<int> stuff2;
        std::mt19937_64 rng(1000);

        int choice1 = 13, choice2 = 19;
        for (size_t i = 0; i < 100; ++i) {
            stuff1.push_back(rng() % choice1);
            stuff2.push_back(rng() % choice2);
            auto& current = collected[std::make_pair(stuff1.back(), stuff2.back())];
            current.push_back(i);
        }

        auto combined = test_combine_factors(stuff1.size(), std::vector<const int*>{stuff1.data(), stuff2.data()});

        // Reference calculation.
        std::vector<int> expected(stuff1.size());
        std::vector<int> factor1, factor2;
        std::vector<size_t> counts;

        size_t counter = 0;
        for (const auto& p : collected) {
            factor1.push_back(p.first.first);
            factor2.push_back(p.first.second);
            counts.push_back(p.second.size());

            for (auto i : p.second) {
                expected[i] = counter;
            }
            ++counter;
        }

        EXPECT_EQ(expected, combined.second);
        EXPECT_EQ(factor1, combined.first.factors[0]);
        EXPECT_EQ(factor2, combined.first.factors[1]);
        EXPECT_EQ(counts, combined.first.counts);
    }
}
