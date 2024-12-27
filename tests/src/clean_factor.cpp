#include <gtest/gtest.h>

#include <random>

#include "scran_aggregate/clean_factor.hpp"
#include "scran_aggregate/combine_factors.hpp"

template<typename Factor_>
std::pair<std::vector<Factor_>, std::vector<int> > test_clean_factor(size_t n, const Factor_* factor) {
    std::vector<int> cleand(n);
    auto levels = scran_aggregate::clean_factor(n, factor, cleand.data());
    return std::make_pair(std::move(levels), std::move(cleand));
}

TEST(CleanFactors, Simple) {
    // Simple sorted case.
    {
        std::vector<int> stuff{0, 0, 1, 1, 1, 2, 2, 2, 2 };
        auto cleand = test_clean_factor(stuff.size(), stuff.data());
        EXPECT_EQ(cleand.second, stuff);

        std::vector<int> expected { 0, 1, 2 };
        EXPECT_EQ(cleand.first, expected);
    }

    // Testing the unsorted case.
    {
        std::vector<int> stuff{ 1, 0, 1, 2, 1, 0, 2, 5, 2 };
        auto cleand = test_clean_factor(stuff.size(), stuff.data());
        std::vector<int> cleaned{ 1, 0, 1, 2, 1, 0, 2, 3, 2 };
        EXPECT_EQ(cleand.second, cleaned);

        std::vector<int> expected { 0, 1, 2, 5 };
        EXPECT_EQ(cleand.first, expected);
    }

    // Non-consecutive still works.
    {
        std::vector<int> stuff{ 1, 3, 5, 7, 9 };
        auto cleand = test_clean_factor(stuff.size(), stuff.data());
        std::vector<int> cleaned { 0, 1, 2, 3, 4 };
        EXPECT_EQ(cleand.second, cleaned);

        std::vector<int> levels { 1, 3, 5, 7, 9 };
        EXPECT_EQ(cleand.first, levels);
    }

    {
        std::vector<int> stuff{ 9, 1, 5, 1, 7, 1, 3 };
        auto cleand = test_clean_factor(stuff.size(), stuff.data());
        std::vector<int> cleaned { 4, 0, 2, 0, 3, 0, 1 };
        EXPECT_EQ(cleand.second, cleaned);

        std::vector<int> levels { 1, 3, 5, 7, 9 };
        EXPECT_EQ(cleand.first, levels);
    }
}
