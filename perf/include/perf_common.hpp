#ifndef PERF_COMMON_HPP
#define PERF_COMMON_HPP

extern "C" {
#include "dft.h"
}

#include <iostream>
#include <limits>
#include <random>

static void prepareData(cfloat32_t* src, int length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    for (int idx = 0; idx < length; ++idx) {
        src[idx].re = dist(gen);
        src[idx].im = dist(gen);
    }
}

static void printStat(uint64_t minCycles, uint64_t totalCycles, int length, int iterations, const std::string& funcName = {}) {
    const uint64_t avgCycles = totalCycles / iterations;
    const float cpeAvg = static_cast<float>(avgCycles) / length;

    const float cpeMin = static_cast<float>(minCycles) / length;

    std::cout << "================================" << std::endl;
    if (funcName.empty())
        std::cout << "Testing DFT" << length << std::endl;
    else
        std::cout << "Testing " << funcName << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << minCycles << " cycles min" << std::endl;
    std::cout << cpeMin << " CPE min" << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << avgCycles << " cycles average" << std::endl;
    std::cout << cpeAvg << " CPE average" << std::endl;
    std::cout << "================================" << std::endl;
}

#ifdef __riscv
static inline uint64_t __attribute__((__gnu_inline__, __always_inline__, __artificial__)) rdcycle(void)
{
    uint64_t dst;
    asm volatile ("csrrs %0, 0xc00, x0" : "=r" (dst));
    return dst;
}
#endif // __riscv

#ifdef __x86_64__
static __inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) uint64_t rdcycle(void){
#ifdef C_MSVC
  return __rdtsc();
#else
  uint64_t a, d;

  __asm__ __volatile__ ("rdtsc" : "=a" (a), "=d" (d));

  return ((uint64_t)a + ((uint64_t)d << 32));
#endif
}
#endif // __x86_64__

#endif // PERF_COMMON_H
