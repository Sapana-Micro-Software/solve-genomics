# Test Suite Documentation

This directory contains comprehensive tests for the DNA sequence alignment algorithms.

## Test Types

### 1. Unit Tests
Test individual functions and methods in isolation:
- `utils_test.cpp` - Utility function tests
- `ExactMatch_test.cpp` - Exact match algorithm tests
- `NaiveSearch_test.cpp` - Naive search algorithm tests
- `SmithWaterman_test.cpp` - Smith-Waterman algorithm tests
- `NeedlemanWunsch_test.cpp` - Needleman-Wunsch algorithm tests

### 2. Regression Tests
Test known good outputs to ensure algorithms haven't regressed:
- Included in each algorithm test file
- Test cases with expected outputs from literature
- Consistency checks across multiple runs

### 3. A-B Comparison Tests
Compare different algorithms on the same inputs:
- `algorithm_comparison_test.cpp` - Compares:
  - ExactMatch vs NaiveSearch (should produce identical results)
  - Smith-Waterman vs Needleman-Wunsch (local vs global alignment)
  - Performance comparisons
  - Different scoring systems

### 4. Integration Tests
Test components working together:
- `integration_test.cpp` - Tests:
  - Complete workflows (read, clean, search, align)
  - File I/O integration
  - Multiple algorithms on same data
  - Error handling across components

### 5. UX Tests
Test user experience and usability:
- Included in `integration_test.cpp`:
  - Consistency and predictability
  - Performance with reasonable input sizes
  - Clear error handling
  - Result interpretability
  - Case-insensitive behavior
  - Formatting and display
  - Edge cases
  - Configurability
  - Real-world scenarios
  - Stability with large inputs

### 6. Blackbox Tests
Test without knowing internal implementation:
- Included in each algorithm test file
- Test properties that should hold regardless of implementation:
  - Output validity
  - Consistency properties
  - Boundary conditions
  - Mathematical properties

## Building and Running Tests

### Prerequisites
- CMake 3.10 or higher
- C++11 compatible compiler
- Google Test (automatically downloaded by CMake)

### Build Tests

```bash
mkdir build
cd build
cmake ..
make
```

### Run All Tests

```bash
cd build
./dna_aligner_tests
```

### Run Specific Test Suites

```bash
# Run only unit tests
./dna_aligner_tests --gtest_filter="*UnitTest*"

# Run only regression tests
./dna_aligner_tests --gtest_filter="*Regression*"

# Run only A-B comparison tests
./dna_aligner_tests --gtest_filter="*Comparison*"

# Run only integration tests
./dna_aligner_tests --gtest_filter="*Integration*"

# Run only UX tests
./dna_aligner_tests --gtest_filter="*UX*"
```

### Run Specific Test Cases

```bash
# Run a specific test
./dna_aligner_tests --gtest_filter="ExactMatchUnitTest.Search_SingleMatch"

# Run all ExactMatch tests
./dna_aligner_tests --gtest_filter="ExactMatch*"

# Run all Smith-Waterman tests
./dna_aligner_tests --gtest_filter="SmithWaterman*"
```

### Using CTest

```bash
cd build
ctest
ctest --verbose
ctest --output-on-failure
```

## Test Coverage

### ExactMatch Tests
- ✅ Single and multiple matches
- ✅ No match cases
- ✅ Empty input handling
- ✅ Case-insensitive matching
- ✅ Pattern longer than sequence
- ✅ Overlapping patterns
- ✅ Consistency checks
- ✅ Position validation

### NaiveSearch Tests
- ✅ All ExactMatch test cases (should match)
- ✅ A-B comparison with ExactMatch
- ✅ Performance characteristics

### Smith-Waterman Tests
- ✅ Identical sequences
- ✅ Similar sequences
- ✅ Sequences with mismatches
- ✅ Empty input handling
- ✅ Custom scoring
- ✅ Position validation
- ✅ Local alignment properties
- ✅ Score non-negativity
- ✅ Aligned sequence properties

### Needleman-Wunsch Tests
- ✅ Identical sequences
- ✅ Sequences with mismatches
- ✅ Sequences with gaps
- ✅ Different length sequences
- ✅ Global alignment properties
- ✅ Empty sequence handling
- ✅ Length conservation
- ✅ Score properties

### Utility Function Tests
- ✅ DNA sequence validation
- ✅ Case conversion
- ✅ Whitespace removal
- ✅ Sequence cleaning
- ✅ File I/O
- ✅ Alignment formatting
- ✅ Identity calculation

## Test Statistics

- **Total Test Files**: 7
- **Total Test Cases**: 100+ individual tests
- **Coverage**: All public APIs and major code paths
- **Test Types**: Unit, Regression, A-B, Integration, UX, Blackbox

## Adding New Tests

When adding new functionality:

1. **Unit Tests**: Add to the appropriate algorithm test file
2. **Regression Tests**: Add known good outputs
3. **A-B Tests**: Compare with existing algorithms if applicable
4. **Integration Tests**: Test with other components
5. **UX Tests**: Ensure user-friendly behavior
6. **Blackbox Tests**: Test expected properties

## Continuous Integration

These tests are designed to run in CI/CD pipelines:

```bash
# Example CI script
mkdir build && cd build
cmake ..
make
./dna_aligner_tests --gtest_output=xml:test_results.xml
```

## Troubleshooting

### Tests Fail to Build
- Ensure CMake can download Google Test (check internet connection)
- Verify C++11 compiler support
- Check that all source files are in the correct directories

### Tests Fail to Run
- Check that the main executable builds successfully first
- Verify all dependencies are linked correctly
- Check for missing include files

### Performance Tests Fail
- Performance tests have timing thresholds that may vary by machine
- Adjust timing expectations if running on slower hardware
- Focus on correctness rather than exact timing

## Test Maintenance

- Update regression tests when algorithm behavior changes intentionally
- Add new test cases when bugs are discovered
- Keep A-B tests synchronized when comparing algorithms
- Review and update UX tests based on user feedback

