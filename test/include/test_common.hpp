#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

extern "C" {
#include "dft.h"
}

#include "gtest/gtest.h"

#include <cmath>
#include <random>
#include <vector>

#define EXPECT_NEAR_VECTOR_VECTOR_COMPLEX(vectorA, vectorB, threshold) \
for (size_t idx{0U}; idx < vectorA.size(); ++idx) \
{ \
    EXPECT_NEAR(vectorA[idx].re, vectorB[idx].re, threshold); \
    EXPECT_NEAR(vectorA[idx].im, vectorB[idx].im, threshold); \
}

[[maybe_unused]]
static void printVector(std::vector<cfloat32_t> vector) {
    for (size_t idx = 0; idx < vector.size(); ++idx) {
        std::cout << vector[idx].re << " " << vector[idx].im << std::endl;
    }
}

#endif // TEST_COMMON_H
