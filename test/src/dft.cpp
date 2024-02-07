#include "test_common.hpp"

#define TEST_2PI (float)0x1.921fb54442d18p+2

class Dft : public ::testing::Test
{
public:
    std::vector<cfloat32_t> source;
    std::vector<cfloat32_t> reference;
    std::vector<cfloat32_t> result;

    void prepareData(int length) {
        source.resize(length);
        reference.resize(length);
        result.resize(length);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        for (int idx = 0; idx < length; ++idx) {
            float phase = TEST_2PI * dist(gen);
            // signal(x, ф) = e^(ix + ф) + e^(3ix + ф) + e^(7ix + ф)
            source[idx].re = cosf(1.0f * TEST_2PI * (float)idx / (float)(4 * length) + phase)
                           + cosf(3.0f * TEST_2PI * (float)idx / (float)(4 * length) + phase)
                           + cosf(7.0f * TEST_2PI * (float)idx / (float)(4 * length) + phase);
            source[idx].im =-sinf(1.0f * TEST_2PI * (float)idx / (float)(4 * length) + phase)
                           + sinf(3.0f * TEST_2PI * (float)idx / (float)(4 * length) + phase)
                           + sinf(7.0f * TEST_2PI * (float)idx / (float)(4 * length) + phase);

            reference[idx].re = 0.0f;
            reference[idx].im = 0.0f;
            result[idx].re = 0.0f;
            result[idx].im = 0.0f;
        }
    }
};

TEST_F(Dft, test_dft2)
{
    static constexpr int length = 2;
    static constexpr float threshold = 1e-2f;

    Dft::prepareData(length);

    refDftFwd(source.data(), reference.data(), length);

    dft2Fwd(source.data(), result.data());

    EXPECT_NEAR_VECTOR_VECTOR_COMPLEX(result, reference, threshold);
}

TEST_F(Dft, test_dft4)
{
    static constexpr int length = 4;
    static constexpr float threshold = 1e-2f;

    Dft::prepareData(length);

    refDftFwd(source.data(), reference.data(), length);

    dft4Fwd(source.data(), result.data());

    EXPECT_NEAR_VECTOR_VECTOR_COMPLEX(result, reference, threshold);
}

TEST_F(Dft, test_dft8)
{
    static constexpr int length = 8;
    static constexpr float threshold = 1e-2f;

    Dft::prepareData(length);

    refDftFwd(source.data(), reference.data(), length);

    dft8Fwd(source.data(), result.data());

    EXPECT_NEAR_VECTOR_VECTOR_COMPLEX(result, reference, threshold);
}

TEST_F(Dft, test_dft16)
{
    static constexpr int length = 16;
    static constexpr float threshold = 1e-1f;

    Dft::prepareData(length);

    refDftFwd(source.data(), reference.data(), length);
    dft16Fwd(source.data(), result.data());

    EXPECT_NEAR_VECTOR_VECTOR_COMPLEX(result, reference, threshold);
}

TEST_F(Dft, test_dft32)
{
    static constexpr int length = 32;
    static constexpr float threshold = 1e-1f;

    Dft::prepareData(length);

    refDftFwd(source.data(), reference.data(), length);
    dft32Fwd(source.data(), result.data());

    EXPECT_NEAR_VECTOR_VECTOR_COMPLEX(result, reference, threshold);
}
