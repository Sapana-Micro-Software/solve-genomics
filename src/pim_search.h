#ifndef PIM_SEARCH_H
#define PIM_SEARCH_H

#include <string>
#include <vector>
#include <memory>

/**
 * Processing-In-Memory (PIM) optimized sequence search
 * Simulates in-memory processing optimizations for faster search
 */
class PIMSearch {
public:
    PIMSearch();
    
    /**
     * Build PIM-optimized index
     * @param sequences Sequences to index
     */
    void buildIndex(const std::vector<std::string>& sequences);
    
    /**
     * Search with PIM optimizations (vectorized operations)
     * @param pattern Pattern to search for
     * @return Positions where pattern found
     */
    std::vector<size_t> search(const std::string& pattern);
    
    /**
     * Batch search (process multiple patterns in parallel)
     * @param patterns Patterns to search
     * @return Results for each pattern
     */
    std::vector<std::vector<size_t>> batchSearch(const std::vector<std::string>& patterns);
    
    /**
     * Vectorized pattern matching (SIMD-like operations)
     * @param sequence Sequence to search in
     * @param pattern Pattern to find
     * @return Positions
     */
    static std::vector<size_t> vectorizedSearch(const std::string& sequence, 
                                               const std::string& pattern);
    
    /**
     * Cache-optimized search (minimize memory access)
     * @param pattern Pattern to search
     * @return Positions
     */
    std::vector<size_t> cacheOptimizedSearch(const std::string& pattern);
    
    /**
     * Parallel search using multiple threads (simulated)
     * @param pattern Pattern to search
     * @param num_threads Number of threads
     * @return Positions
     */
    std::vector<size_t> parallelSearch(const std::string& pattern, int num_threads = 4);
    
private:
    std::vector<std::string> indexed_sequences_;
    
    /**
     * Precompute pattern hash for fast lookup
     */
    size_t hashPattern(const std::string& pattern);
    
    /**
     * Vectorized comparison (process multiple characters at once)
     */
    bool vectorizedCompare(const std::string& str1, size_t pos1,
                         const std::string& str2, size_t pos2,
                         size_t length);
    
    /**
     * Build suffix array for fast search
     */
    void buildSuffixArray();
    
    std::vector<size_t> suffix_array_;
    std::vector<std::string> suffixes_;
};

#endif // PIM_SEARCH_H

