file(GLOB SOURCES "*/*.cpp")

add_executable(fft_test ${SOURCES})

target_include_directories(fft_test
    PUBLIC
    ${PROJECT_SOURCE_DIR}/test/include
    ${PROJECT_SOURCE_DIR}/3rd-party/gtest/googletest/include
    ${PROJECT_SOURCE_DIR}/lib/include
)

target_link_libraries(fft_test
PUBLIC
    GTest::gtest
    GTest::gtest_main
    dft
PRIVATE
    CommonConfiguration
)


