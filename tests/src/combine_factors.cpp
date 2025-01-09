#include <gtest/gtest.h>

#include <random>

#include "utils.h" // must be included before scran_aggregate
#include "scran_aggregate/combine_factors.hpp"

template<typename Factor_>
std::pair<std::vector<std::vector<Factor_> >, std::vector<int> > test_combine_factors(size_t n, const std::vector<const Factor_*>& factors) {
    std::vector<int> combined(n
#ifdef SCRAN_AGGREGATE_TEST_INIT
        , SCRAN_AGGREGATE_TEST_INIT
#endif
    );
    auto levels = scran_aggregate::combine_factors(n, factors, combined.data());
    return std::make_pair(std::move(levels), std::move(combined));
}

TEST(CombineFactors, Special) {
    // Single factor, just defers to clean_factor.
    {
        std::vector<int> stuff{0, 0, 1, 1, 1, 2, 2, 2, 2 };
        auto combined = test_combine_factors(stuff.size(), std::vector<const int*>{stuff.data()});
        EXPECT_EQ(combined.second, stuff);

        EXPECT_EQ(combined.first.size(), 1);
        std::vector<int> expected { 0, 1, 2 };
        EXPECT_EQ(combined.first[0], expected);
    }

    // Nothing at all.
    {
        auto combined = test_combine_factors(10, std::vector<const int*>{});
        EXPECT_EQ(combined.second, std::vector<int>(10));
        EXPECT_TRUE(combined.first.empty());
    }
}

TEST(CombineFactors, Multiple) {
    // Single factor, the other is a no-op.
    {
        std::vector<int> stuff1{0, 0, 1, 1, 1, 2, 2, 2, 2 };
        std::vector<int> stuff2(stuff1.size());
        auto combined = test_combine_factors(stuff1.size(), std::vector<const int*>{stuff1.data(), stuff2.data()});
        EXPECT_EQ(combined.second, stuff1);

        EXPECT_EQ(combined.first.size(), 2);
        std::vector<int> expected1 { 0, 1, 2 };
        EXPECT_EQ(combined.first[0], expected1);
        std::vector<int> expected2 { 0, 0, 0 };
        EXPECT_EQ(combined.first[1], expected2);
    }

    {
        std::vector<int> stuff1{ 1, 0, 1, 2, 1, 0, 2, 3, 2 };
        std::vector<int> stuff2(stuff1.size());
        auto combined = test_combine_factors(stuff1.size(), std::vector<const int*>{stuff1.data(), stuff2.data()});
        EXPECT_EQ(combined.second, stuff1);

        EXPECT_EQ(combined.first.size(), 2);
        std::vector<int> expected1 { 0, 1, 2, 3 };
        EXPECT_EQ(combined.first[0], expected1);
        std::vector<int> expected2 { 0, 0, 0, 0 };
        EXPECT_EQ(combined.first[1], expected2);
    }

    {
        std::vector<int> stuff2{ 1, 3, 5, 7, 9 };
        std::vector<int> stuff1(stuff2.size());
        auto combined = test_combine_factors(stuff1.size(), std::vector<const int*>{stuff1.data(), stuff2.data()});
        std::vector<int> expected { 0, 1, 2, 3, 4 };
        EXPECT_EQ(combined.second, expected);

        EXPECT_EQ(combined.first.size(), 2);
        std::vector<int> expected1 { 0, 0, 0, 0, 0 };
        EXPECT_EQ(combined.first[0], expected1);
        std::vector<int> expected2 { 1, 3, 5, 7, 9 };
        EXPECT_EQ(combined.first[1], expected2);
    }

    // Multiple non-trivial factors.
    {
        std::vector<int> stuff1{ 0, 0, 1, 1, 1, 2, 2, 2, 2 };
        std::vector<int> stuff2{ 0, 1, 2, 0, 1, 2, 0, 1, 2 };
        auto combined = test_combine_factors(stuff1.size(), std::vector<const int*>{stuff1.data(), stuff2.data()});

        std::vector<int> expected { 0, 1, 4, 2, 3, 7, 5, 6, 7 };
        EXPECT_EQ(combined.second, expected);

        EXPECT_EQ(combined.first.size(), 2);
        std::vector<int> levels1 { 0, 0, 1, 1, 1, 2, 2, 2 };
        std::vector<int> levels2 { 0, 1, 0, 1, 2, 0, 1, 2 };
        EXPECT_EQ(combined.first[0], levels1);
        EXPECT_EQ(combined.first[1], levels2);
    }

    // Simulated example.
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
        EXPECT_EQ(factor1, combined.first[0]);
        EXPECT_EQ(factor2, combined.first[1]);
    }
}

template<typename Factor_, typename Number_>
std::pair<std::vector<std::vector<Factor_> >, std::vector<int> > test_combine_factors_unused(size_t n, const std::vector<std::pair<const Factor_*, Number_> >& factors) {
    std::vector<int> combined(n
#ifdef SCRAN_AGGREGATE_TEST_INIT
        , SCRAN_AGGREGATE_TEST_INIT
#endif
    );
    auto levels = scran_aggregate::combine_factors_unused(n, factors, combined.data());
    return std::make_pair(std::move(levels), std::move(combined));
}

TEST(CombineFactorsUnused, Special) {
    std::vector<int> stuff{ 1, 3, 5, 3, 1 };
    auto combined = test_combine_factors_unused(stuff.size(), std::vector<std::pair<const int*, int> >{ { stuff.data(), 7 } });
    EXPECT_EQ(combined.second, stuff);
    std::vector<int> levels{ 0, 1, 2, 3, 4, 5, 6 };
    EXPECT_EQ(combined.first[0], levels);

    auto combined2 = test_combine_factors_unused(10, std::vector<std::pair<const int*, int> >{});
    EXPECT_EQ(combined2.second, std::vector<int>(10));
    EXPECT_TRUE(combined2.first.empty());
}

TEST(CombineFactorsUnused, Multiple) {
    std::vector<int> stuff1{ 0, 0, 1, 1, 1, 2, 2, 2, 2 };
    std::vector<int> stuff2{ 0, 1, 2, 0, 1, 2, 0, 1, 2 };

    {
        auto combined = test_combine_factors_unused(stuff1.size(), 
            std::vector<std::pair<const int*, int> >{ 
                { stuff1.data(), 3 },
                { stuff2.data(), 3 }
            }
        );

        std::vector<int> expected { 0, 1, 5, 3, 4, 8, 6, 7, 8 };
        EXPECT_EQ(combined.second, expected);

        EXPECT_EQ(combined.first.size(), 2);
        std::vector<int> levels1 { 0, 0, 0, 1, 1, 1, 2, 2, 2 };
        std::vector<int> levels2 { 0, 1, 2, 0, 1, 2, 0, 1, 2 }; // (0, 2) is included even though it is not observed.
        EXPECT_EQ(combined.first[0], levels1);
        EXPECT_EQ(combined.first[1], levels2);
    }

    {
        auto combined = test_combine_factors_unused(stuff1.size(), 
            std::vector<std::pair<const int*, int> >{ 
                { stuff1.data(), 4 },
                { stuff2.data(), 3 }
            }
        );

        std::vector<int> expected { 0, 1, 5, 3, 4, 8, 6, 7, 8 };
        EXPECT_EQ(combined.second, expected);

        EXPECT_EQ(combined.first.size(), 2);
        std::vector<int> levels1 { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 };
        std::vector<int> levels2 { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 }; // (0, 2) and (3, *) are included even though they are not observed.
        EXPECT_EQ(combined.first[0], levels1);
        EXPECT_EQ(combined.first[1], levels2);
    }

    // Multiple things at play here.
    {
        std::vector<int> mock{ 1 };
        auto combined = test_combine_factors_unused(1,
            std::vector<std::pair<const int*, int> >{ 
                { mock.data(), 2 },
                { mock.data(), 3 },
                { mock.data(), 4 }
            }
        );

        EXPECT_EQ(combined.second[0], 17); // i.e., 1*3*4 + 1*4 + 1.

        auto create_mock_sequence = [](size_t nlevels, size_t inner, size_t outer) -> auto { 
            std::vector<int> output(nlevels * inner * outer);
            auto oIt = output.begin();
            for (size_t o = 0; o < outer; ++o) {
                for (size_t l = 0; l < nlevels; ++l) {
                    std::fill_n(oIt, inner, l);
                    oIt += inner;
                }
            }
            return output;
        };

        EXPECT_EQ(combined.first[0], create_mock_sequence(2, 12, 1));
        EXPECT_EQ(combined.first[1], create_mock_sequence(3, 4, 2));
        EXPECT_EQ(combined.first[2], create_mock_sequence(4, 1, 6));
    }
}
