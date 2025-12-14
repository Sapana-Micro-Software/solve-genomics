# Benchmark Suite Documentation

## Overview

The benchmark suite tests DNA sequence alignment algorithms on sequences with different complexity characteristics:
- **High Entropy**: Random sequences with maximum information content
- **Low Entropy**: Highly repetitive sequences with minimal information content
- **Moderate Complexity**: Mixed sequences with both random and repetitive regions

## Building Benchmarks

```bash
mkdir build
cd build
cmake ..
make
```

This creates the `dna_benchmark` executable in the `bin/` directory.

## Running Benchmarks

### Standard Benchmark Suite

```bash
./bin/dna_benchmark
```

This runs comprehensive benchmarks with:
- Multiple sequence lengths: 100, 500, 1000, 5000, 10000 bases
- Multiple pattern lengths: 4, 8, 16, 32 bases
- All search algorithms: ExactMatch, NaiveSearch, FuzzySearch (with different distance thresholds)
- Alignment algorithms: Smith-Waterman, Needleman-Wunsch (for sequences ≤ 1000 bases)

### Detailed Benchmark

```bash
./bin/dna_benchmark detailed
```

This runs a detailed benchmark with:
- Fixed sequence length (1000 bases)
- Multiple iterations (10) for statistical significance
- Focused analysis on different sequence types

## Sequence Types

### High Entropy Sequences
- **Characteristics**: Random distribution of A, T, G, C
- **Entropy**: ~2.0 bits (maximum for 4-letter alphabet)
- **Use Case**: Testing worst-case performance, random data
- **Example**: `ATCGGCATAGCTAGCTAG...` (random)

### Low Entropy Sequences
- **Characteristics**: Highly repetitive patterns
- **Entropy**: < 1.0 bits (low information content)
- **Use Case**: Testing best-case performance, repetitive regions
- **Example**: `ATCGATCGATCGATCG...` (repeating pattern)

### Tandem Repeats
- **Characteristics**: Specific unit repeated many times
- **Use Case**: Testing on common genomic feature
- **Example**: `ATCGATCGATCG...` (tandem repeat of "ATCG")

### Moderate Complexity
- **Characteristics**: Mix of random and repetitive regions
- **Use Case**: Testing on realistic genomic sequences
- **Example**: Combination of repetitive and random segments

## Metrics Collected

### Performance Metrics
- **Execution Time**: Measured in milliseconds (ms)
- **Algorithm**: Algorithm name
- **Sequence Type**: HighEntropy, LowEntropy, Moderate, TandemRepeat

### Sequence Characteristics
- **Sequence Length**: Length of the sequence being searched
- **Pattern Length**: Length of the search pattern
- **Entropy**: Shannon entropy in bits (measure of randomness)
- **GC Content**: Ratio of G and C bases (0.0 to 1.0)
- **Matches Found**: Number of matches/alignments found

## Output Formats

### Console Output
Results are displayed in a formatted table showing:
- Algorithm name
- Sequence type
- Sequence and pattern lengths
- Execution time
- Number of matches
- Entropy and GC content

### CSV Output
Results are automatically saved to `benchmark_results.csv` with columns:
- Algorithm
- SequenceType
- SequenceLength
- PatternLength
- Time_ms
- Matches
- Entropy
- GCContent

## Expected Results

### Performance Characteristics

1. **ExactMatch vs NaiveSearch**
   - Should have similar performance (both O(n*m))
   - ExactMatch may be slightly faster due to optimizations

2. **FuzzySearch**
   - Performance decreases with increasing distance threshold
   - FuzzySearch_0 (distance 0) ≈ ExactMatch performance
   - Higher thresholds (1, 2) are slower due to edit distance calculation

3. **Sequence Complexity Impact**
   - **High Entropy**: Generally slower (more comparisons needed)
   - **Low Entropy**: Faster (early matches in repetitive regions)
   - **Moderate**: Performance between high and low entropy

4. **Alignment Algorithms**
   - Smith-Waterman and Needleman-Wunsch: O(n*m) complexity
   - Performance scales quadratically with sequence length
   - Low entropy may show different alignment patterns

### Typical Performance Ranges

| Algorithm | Sequence Length | Pattern Length | Time (ms) |
|-----------|----------------|----------------|-----------|
| ExactMatch | 1000 | 10 | 0.1 - 1.0 |
| NaiveSearch | 1000 | 10 | 0.1 - 1.0 |
| FuzzySearch_1 | 1000 | 10 | 1.0 - 10.0 |
| SmithWaterman | 1000 | 1000 | 10 - 100 |
| NeedlemanWunsch | 1000 | 1000 | 10 - 100 |

*Note: Actual times vary based on hardware and sequence characteristics*

## Analyzing Results

### Using the CSV File

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv('benchmark_results.csv')

# Compare algorithms
algorithms = df.groupby('Algorithm')['Time_ms'].mean()
print(algorithms)

# Compare by sequence type
by_type = df.groupby(['Algorithm', 'SequenceType'])['Time_ms'].mean()
print(by_type)

# Plot performance vs sequence length
for algo in df['Algorithm'].unique():
    subset = df[df['Algorithm'] == algo]
    plt.plot(subset['SequenceLength'], subset['Time_ms'], label=algo)
plt.xlabel('Sequence Length')
plt.ylabel('Time (ms)')
plt.legend()
plt.show()
```

### Key Insights

1. **Algorithm Selection**: Use results to choose the right algorithm for your use case
2. **Performance Scaling**: Understand how algorithms scale with input size
3. **Complexity Impact**: See how sequence complexity affects performance
4. **Trade-offs**: Balance between accuracy (fuzzy search) and speed (exact match)

## Custom Benchmarks

You can modify `src/benchmark.cpp` to:
- Test different sequence lengths
- Test custom patterns
- Add new sequence types
- Measure additional metrics (memory usage, cache performance)
- Compare specific algorithms

## Reproducibility

The benchmark uses a fixed random seed (42) for reproducibility. To change:
- Modify the seed in `BenchmarkRunner` constructor
- Or use `SequenceGenerator` with a custom seed

## Troubleshooting

### Benchmarks Run Slowly
- Reduce sequence lengths in `runComplexityBenchmarks()`
- Skip alignment algorithms for very long sequences
- Run detailed benchmark instead of full suite

### Out of Memory
- Reduce maximum sequence length
- Skip alignment algorithms (they use O(n*m) memory)

### Inconsistent Results
- Run multiple iterations (detailed benchmark does this)
- Check for system load (other processes)
- Ensure sufficient system resources

## Future Enhancements

Potential additions to the benchmark suite:
- Memory usage tracking
- Cache performance analysis
- Multi-threaded algorithm comparisons
- Real genomic sequence datasets
- Statistical analysis (mean, std dev, confidence intervals)
- Visualization tools
- Performance regression testing

