#ifndef CONCURRENT_MULTI_SEARCH_H
#define CONCURRENT_MULTI_SEARCH_H

#include "SequenceAligner.h"
#include "ExactMatch.h"
#include "NaiveSearch.h"
#include "FuzzySearch.h"
#include "classic_string_matching.h"
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <future>
#include <atomic>
#include <memory>

/**
 * Concurrent Multi-Technique Substring Search
 * Runs multiple search algorithms simultaneously using separate threads
 * and combines their results for comprehensive pattern matching
 */
class ConcurrentMultiSearch {
public:
    /**
     * Search technique enumeration
     */
    enum SearchTechnique {
        EXACT_MATCH,           // ExactMatch algorithm
        NAIVE_SEARCH,          // NaiveSearch algorithm
        RABIN_KARP,            // Rabin-Karp rolling hash
        KMP,                   // Knuth-Morris-Pratt
        BOYER_MOORE,           // Boyer-Moore heuristic
        FUZZY_SEARCH,          // Fuzzy search with edit distance
        FUZZY_HAMMING,         // Fuzzy search with Hamming distance
        ALL_TECHNIQUES         // Use all available techniques
    };
    
    /**
     * Combined search result from multiple techniques
     */
    struct MultiTechniqueResult {
        std::map<std::string, std::vector<int>> technique_positions;  // Positions by technique name
        std::map<std::string, int> technique_counts;                 // Match counts by technique
        std::map<std::string, double> technique_times;                // Execution times (microseconds)
        std::vector<int> consensus_positions;                         // Positions found by multiple techniques
        std::vector<int> unique_positions;                             // Positions found by only one technique
        int total_unique_matches;                                      // Total unique match positions
        double total_time;                                             // Total execution time
        int num_techniques_used;                                       // Number of techniques executed
        
        MultiTechniqueResult() : total_unique_matches(0), total_time(0.0), num_techniques_used(0) {}
    };
    
    /**
     * Configuration for concurrent search
     */
    struct SearchConfig {
        std::vector<SearchTechnique> techniques;                       // Techniques to use
        int max_errors;                                                // Maximum errors for fuzzy search
        bool require_consensus;                                        // Require multiple techniques to agree
        int min_consensus_count;                                       // Minimum techniques that must agree
        bool use_threads;                                              // Use threads (true) or processes (false)
        int num_threads;                                              // Number of threads (0 = auto)
        
        SearchConfig() : max_errors(0), require_consensus(false), 
                        min_consensus_count(2), use_threads(true), num_threads(0) {}
    };
    
    ConcurrentMultiSearch();
    
    /**
     * Search using multiple techniques concurrently
     * @param sequence DNA sequence to search in
     * @param pattern Pattern to find
     * @param config Search configuration
     * @return MultiTechniqueResult with combined results
     */
    MultiTechniqueResult search(const std::string& sequence,
                               const std::string& pattern,
                               const SearchConfig& config = SearchConfig());
    
    /**
     * Search using specified techniques
     * @param sequence DNA sequence to search in
     * @param pattern Pattern to find
     * @param techniques Vector of techniques to use
     * @return MultiTechniqueResult with combined results
     */
    MultiTechniqueResult search(const std::string& sequence,
                               const std::string& pattern,
                               const std::vector<SearchTechnique>& techniques);
    
    /**
     * Get default configuration (all techniques)
     */
    static SearchConfig getDefaultConfig();
    
    /**
     * Get fast configuration (exact techniques only)
     */
    static SearchConfig getFastConfig();
    
    /**
     * Get comprehensive configuration (all techniques including fuzzy)
     */
    static SearchConfig getComprehensiveConfig();
    
private:
    /**
     * Search result from a single technique
     */
    struct TechniqueResult {
        std::string technique_name;                                   // Name of technique
        std::vector<int> positions;                                    // Found positions
        int count;                                                     // Number of matches
        double execution_time;                                         // Execution time in microseconds
        bool success;                                                  // Whether search succeeded
        
        TechniqueResult() : count(0), execution_time(0.0), success(false) {}
    };
    
    /**
     * Execute a single search technique
     */
    TechniqueResult executeTechnique(SearchTechnique technique,
                                    const std::string& sequence,
                                    const std::string& pattern,
                                    int max_errors);
    
    /**
     * Execute search using threads
     */
    MultiTechniqueResult searchWithThreads(const std::string& sequence,
                                           const std::string& pattern,
                                           const SearchConfig& config);
    
    /**
     * Execute search using processes (future implementation)
     */
    MultiTechniqueResult searchWithProcesses(const std::string& sequence,
                                            const std::string& pattern,
                                            const SearchConfig& config);
    
    /**
     * Combine results from multiple techniques
     */
    void combineResults(const std::vector<TechniqueResult>& results,
                       MultiTechniqueResult& combined_result,
                       const SearchConfig& config);
    
    /**
     * Find consensus positions (found by multiple techniques)
     */
    std::vector<int> findConsensusPositions(
            const std::map<std::string, std::vector<int>>& technique_positions,
            int min_consensus);
    
    /**
     * Find unique positions (found by only one technique)
     */
    std::vector<int> findUniquePositions(
            const std::map<std::string, std::vector<int>>& technique_positions);
    
    /**
     * Get technique name
     */
    std::string getTechniqueName(SearchTechnique technique);
    
    /**
     * Mutex for thread-safe result collection
     */
    std::mutex result_mutex_;
    
    /**
     * Atomic counter for completed searches
     */
    std::atomic<int> completed_searches_;
};

#endif // CONCURRENT_MULTI_SEARCH_H

