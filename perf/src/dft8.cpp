#include "perf_common.hpp"

int main() {
    constexpr int length = 8;
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
        dft8Fwd(src, dst);
        end = rdcycle();
        diff = end - start;
        totalCycles += diff;
        if (diff < minCycles)
            minCycles = diff;
    }

    printStat(minCycles, totalCycles, length, iterations);

    return 0;
}
