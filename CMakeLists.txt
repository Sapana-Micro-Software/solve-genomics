cmake_minimum_required(VERSION 3.10)
project(DNAAlignment VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# Source files
set(SOURCES
    src/main.cpp
    src/utils.cpp
    src/ExactMatch.cpp
    src/NaiveSearch.cpp
    src/SmithWaterman.cpp
    src/NeedlemanWunsch.cpp
    src/FuzzySearch.cpp
    src/grammar_compression.cpp
    src/association_list.cpp
    src/lossy_compression.cpp
    src/edit_distance_algorithms.cpp
    src/embedding_search.cpp
    src/cnn_sequence.cpp
    src/lightweight_llm.cpp
    src/pim_search.cpp
    src/classic_string_matching.cpp
    src/mcmc_search.cpp
    src/ddmcmc_search.cpp
    src/parallel_search.cpp
    src/warp_ctc.cpp
    src/aho_corasick.cpp
    src/suffix_tree.cpp
    src/suffix_array.cpp
    src/advanced_dp.cpp
    src/deep_learning_advanced.cpp
    src/wu_manber.cpp
    src/bitap.cpp
    src/graph_based.cpp
    src/heuristic_search.cpp
    src/vector_similarity.cpp
    src/variant_calling.cpp
    src/metagenomics.cpp
    src/chord_hashing.cpp
    src/sequence_generator.cpp
    src/probabilistic_ml.cpp
    src/concurrent_multi_search.cpp
    src/skip_graph.cpp
    src/dancing_links.cpp
)

# Header files for sequence generator (used by benchmark)
set(SEQUENCE_GEN_HEADERS
    src/sequence_generator.h
)

# Header files
set(HEADERS
    src/SequenceAligner.h
    src/utils.h
    src/ExactMatch.h
    src/NaiveSearch.h
    src/SmithWaterman.h
    src/NeedlemanWunsch.h
    src/FuzzySearch.h
    src/compression.h
    src/grammar_compression.h
    src/association_list.h
    src/lossy_compression.h
    src/edit_distance_algorithms.h
    src/embedding_search.h
    src/cnn_sequence.h
    src/lightweight_llm.h
    src/pim_search.h
    src/classic_string_matching.h
    src/mcmc_search.h
    src/ddmcmc_search.h
    src/parallel_search.h
    src/warp_ctc.h
    src/aho_corasick.h
    src/suffix_tree.h
    src/suffix_array.h
    src/advanced_dp.h
    src/deep_learning_advanced.h
    src/wu_manber.h
    src/bitap.h
    src/graph_based.h
    src/heuristic_search.h
    src/vector_similarity.h
    src/variant_calling.h
    src/metagenomics.h
    src/chord_hashing.h
    src/probabilistic_ml.h
    src/concurrent_multi_search.h
    src/skip_graph.h
    src/dancing_links.h
)

# Create executable
add_executable(dna_aligner ${SOURCES} ${HEADERS})

# Set output directory
set_target_properties(dna_aligner PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Benchmark executable
set(BENCHMARK_SOURCES
    src/benchmark.cpp
    src/sequence_generator.cpp
    src/utils.cpp
    src/ExactMatch.cpp
    src/NaiveSearch.cpp
    src/FuzzySearch.cpp
    src/SmithWaterman.cpp
    src/NeedlemanWunsch.cpp
    src/grammar_compression.cpp
    src/association_list.cpp
    src/lossy_compression.cpp
)

add_executable(dna_benchmark ${BENCHMARK_SOURCES} ${HEADERS} src/sequence_generator.h)

set_target_properties(dna_benchmark PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Installation (optional)
install(TARGETS dna_aligner
    RUNTIME DESTINATION bin
)

# Print configuration
message(STATUS "DNA Alignment Tool Configuration:")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# Testing configuration
enable_testing()

# Option to use system Google Test or download it
option(FETCH_GOOGLETEST "Download Google Test from GitHub" ON)

if(FETCH_GOOGLETEST)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
else()
    find_package(GTest REQUIRED)
endif()

# Test source files
set(TEST_SOURCES
    tests/test_main.cpp
    tests/utils_test.cpp
    tests/ExactMatch_test.cpp
    tests/NaiveSearch_test.cpp
    tests/SmithWaterman_test.cpp
    tests/NeedlemanWunsch_test.cpp
    tests/FuzzySearch_test.cpp
    tests/algorithm_comparison_test.cpp
    tests/integration_test.cpp
)

# Test library sources (without main.cpp)
set(TEST_LIB_SOURCES
    src/utils.cpp
    src/ExactMatch.cpp
    src/NaiveSearch.cpp
    src/SmithWaterman.cpp
    src/NeedlemanWunsch.cpp
    src/FuzzySearch.cpp
    src/grammar_compression.cpp
    src/association_list.cpp
    src/lossy_compression.cpp
    src/sequence_generator.cpp
)

# Create test executable
add_executable(dna_aligner_tests ${TEST_SOURCES} ${TEST_LIB_SOURCES})

# Link test executable with Google Test
target_link_libraries(dna_aligner_tests
    PRIVATE
    GTest::gtest
    GTest::gtest_main
    GTest::gmock
)

# Include directories for tests
target_include_directories(dna_aligner_tests PRIVATE src)

# Add test to CTest
add_test(NAME AllTests COMMAND dna_aligner_tests)

