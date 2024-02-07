#include "perf_common.hpp"

int main() {
    constexpr int length = 6;
    constexpr int iterations = 100'000;

    volatile uint64_t start = 0,
                      end = 0,
                      totalCycles = 0,
                      minCycles = std::numeric_limits<uint64_t>::max(),
                      diff;

    cfloat32_t src[length];
    cfloat32_t dst[length];

    prepareData(src, length);

    for(int i = 0; i < iterations; i++){
        start = rdcycle();
        refDftFwd(src, dst, length);
        end = rdcycle();
        diff = end - start;
        totalCycles += diff;
        if (diff < minCycles)
            minCycles = diff;
    }

    std::string funcName = std::string{"reference DFT with length "} + std::to_string(length);

    printStat(minCycles, totalCycles, length, iterations, funcName);

    return 0;
}
