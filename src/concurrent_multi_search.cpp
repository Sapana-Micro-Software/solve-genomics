#include "concurrent_multi_search.h"
#include <algorithm>
#include <chrono>
#include <thread>
#include <future>
#include <set>
#include <sstream>

ConcurrentMultiSearch::ConcurrentMultiSearch() : completed_searches_(0) {
}

ConcurrentMultiSearch::MultiTechniqueResult ConcurrentMultiSearch::search(
        const std::string& sequence,
        const std::string& pattern,
        const SearchConfig& config) {
    if (sequence.empty() || pattern.empty()) {                                                  // Validate input parameters: check if sequences are non-empty
        return MultiTechniqueResult();                                                         // Return empty result structure if invalid input
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();                              // Record start time for total execution measurement
    
    MultiTechniqueResult result;                                                                // Initialize combined result structure
    
    if (config.use_threads) {                                                                  // Check if configuration specifies thread-based execution
        result = searchWithThreads(sequence, pattern, config);                                 // Execute search using multiple threads
    } else {
        result = searchWithProcesses(sequence, pattern, config);                               // Execute search using multiple processes (future implementation)
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();                                // Record end time
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); // Calculate total execution duration
    result.total_time = duration.count();                                                      // Store total execution time in microseconds
    
    return result;                                                                              // Return combined result from all techniques
}

ConcurrentMultiSearch::MultiTechniqueResult ConcurrentMultiSearch::search(
        const std::string& sequence,
        const std::string& pattern,
        const std::vector<SearchTechnique>& techniques) {
    SearchConfig config;                                                                       // Create default search configuration
    config.techniques = techniques;                                                             // Set list of techniques to use
    return search(sequence, pattern, config);                                                  // Call main search method with configuration
}

ConcurrentMultiSearch::SearchConfig ConcurrentMultiSearch::getDefaultConfig() {
    SearchConfig config;                                                                       // Create default search configuration
    config.techniques = {EXACT_MATCH, NAIVE_SEARCH, KMP, BOYER_MOORE};                       // Use fast exact matching techniques
    config.use_threads = true;                                                                 // Enable thread-based concurrent execution
    config.num_threads = 0;                                                                    // Auto-detect optimal thread count from hardware
    config.require_consensus = false;                                                          // Don't require multiple techniques to agree
    return config;                                                                              // Return default configuration
}

ConcurrentMultiSearch::SearchConfig ConcurrentMultiSearch::getFastConfig() {
    SearchConfig config;                                                                       // Create fast search configuration
    config.techniques = {KMP, BOYER_MOORE};                                                  // Use only the fastest exact matching techniques
    config.use_threads = true;                                                                 // Enable thread-based concurrent execution
    config.num_threads = 2;                                                                    // Use exactly 2 threads for parallel execution
    config.require_consensus = false;                                                          // Don't require multiple techniques to agree
    return config;                                                                              // Return fast configuration
}

ConcurrentMultiSearch::SearchConfig ConcurrentMultiSearch::getComprehensiveConfig() {
    SearchConfig config;                                                                       // Create comprehensive search configuration
    config.techniques = {EXACT_MATCH, NAIVE_SEARCH, RABIN_KARP, KMP, BOYER_MOORE,            // Use all available search techniques
                        FUZZY_SEARCH, FUZZY_HAMMING};
    config.use_threads = true;                                                                 // Enable thread-based concurrent execution
    config.num_threads = 0;                                                                    // Auto-detect optimal thread count from hardware
    config.require_consensus = true;                                                          // Require multiple techniques to agree on matches
    config.min_consensus_count = 2;                                                            // At least 2 techniques must find same position
    config.max_errors = 1;                                                                     // Allow maximum 1 error for fuzzy search algorithms
    return config;                                                                              // Return comprehensive configuration
}

ConcurrentMultiSearch::MultiTechniqueResult ConcurrentMultiSearch::searchWithThreads(
        const std::string& sequence,
        const std::string& pattern,
        const SearchConfig& config) {
    MultiTechniqueResult result;                                                                // Initialize combined result structure
    
    std::vector<SearchTechnique> techniques = config.techniques;                              // Get list of techniques to use from configuration
    if (techniques.empty() || techniques[0] == ALL_TECHNIQUES) {                               // Check if all techniques requested or list is empty
        techniques = {EXACT_MATCH, NAIVE_SEARCH, RABIN_KARP, KMP, BOYER_MOORE,               // Use all available search techniques
                     FUZZY_SEARCH, FUZZY_HAMMING};
    }
    
    int num_threads = config.num_threads;                                                      // Get configured number of threads
    if (num_threads == 0) {                                                                    // If auto-detect thread count
        num_threads = std::min(static_cast<int>(techniques.size()),                           // Use minimum of: number of techniques
                              static_cast<int>(std::thread::hardware_concurrency()));         // and hardware thread concurrency
    }
    
    std::vector<std::future<TechniqueResult>> futures;                                         // Store future objects for asynchronous execution
    std::vector<TechniqueResult> results;                                                      // Store results from all search techniques
    
    // Launch searches in parallel using async
    for (SearchTechnique technique : techniques) {                                             // Iterate through each search technique
        futures.push_back(std::async(std::launch::async,                                      // Launch asynchronous task (new thread)
                                    [this, technique, &sequence, &pattern, &config]() {       // Lambda function capturing context
            return executeTechnique(technique, sequence, pattern, config.max_errors);         // Execute specific search technique
        }));
    }
    
    // Collect results from all futures
    for (auto& future : futures) {                                                            // Iterate through all future objects
        TechniqueResult tech_result = future.get();                                           // Wait for completion and retrieve result
        results.push_back(tech_result);                                                        // Store technique result
    }
    
    // Combine results
    combineResults(results, result, config);                                                  // Merge results from all techniques into combined result
    
    result.num_techniques_used = static_cast<int>(techniques.size());                         // Store total number of techniques executed
    
    return result;                                                                              // Return combined result from all techniques
}

ConcurrentMultiSearch::MultiTechniqueResult ConcurrentMultiSearch::searchWithProcesses(
        const std::string& sequence,
        const std::string& pattern,
        const SearchConfig& config) {
    // For now, processes are implemented using threads
    // In a full implementation, would use fork() or process pools
    return searchWithThreads(sequence, pattern, config);                                     // Fall back to thread-based implementation
}

ConcurrentMultiSearch::TechniqueResult ConcurrentMultiSearch::executeTechnique(
        SearchTechnique technique,
        const std::string& sequence,
        const std::string& pattern,
        int max_errors) {
    TechniqueResult result;                                                                     // Initialize technique-specific result structure
    result.technique_name = getTechniqueName(technique);                                       // Set human-readable technique name
    
    auto start_time = std::chrono::high_resolution_clock::now();                              // Record start time for this technique
    
    try {
        switch (technique) {                                                                   // Switch on technique enumeration value
            case EXACT_MATCH: {                                                                // Exact match technique
                ExactMatch matcher;                                                            // Create ExactMatch algorithm instance
                SearchResult search_result = matcher.search(sequence, pattern);                // Execute exact substring search
                result.positions = search_result.positions;                                   // Store all match positions found
                result.count = search_result.count;                                            // Store total number of matches
                result.success = true;                                                         // Mark execution as successful
                break;
            }
            case NAIVE_SEARCH: {                                                              // Naive search technique
                NaiveSearch searcher;                                                         // Create NaiveSearch algorithm instance
                SearchResult search_result = searcher.search(sequence, pattern);               // Execute naive brute-force search
                result.positions = search_result.positions;                                   // Store all match positions found
                result.count = search_result.count;                                            // Store total number of matches
                result.success = true;                                                         // Mark execution as successful
                break;
            }
            case RABIN_KARP: {                                                               // Rabin-Karp technique
                ClassicStringMatching::RabinKarp rk;                                         // Create Rabin-Karp rolling hash matcher
                SearchResult search_result = rk.search(sequence, pattern);                    // Execute Rabin-Karp search
                result.positions = search_result.positions;                                   // Store all match positions found
                result.count = search_result.count;                                            // Store total number of matches
                result.success = true;                                                         // Mark execution as successful
                break;
            }
            case KMP: {                                                                      // KMP technique
                ClassicStringMatching::KMP kmp;                                              // Create Knuth-Morris-Pratt matcher
                SearchResult search_result = kmp.search(sequence, pattern);                    // Execute KMP search with failure function
                result.positions = search_result.positions;                                   // Store all match positions found
                result.count = search_result.count;                                            // Store total number of matches
                result.success = true;                                                         // Mark execution as successful
                break;
            }
            case BOYER_MOORE: {                                                              // Boyer-Moore technique
                ClassicStringMatching::BoyerMoore bm;                                         // Create Boyer-Moore heuristic matcher
                SearchResult search_result = bm.search(sequence, pattern);                     // Execute Boyer-Moore search
                result.positions = search_result.positions;                                   // Store all match positions found
                result.count = search_result.count;                                            // Store total number of matches
                result.success = true;                                                         // Mark execution as successful
                break;
            }
            case FUZZY_SEARCH: {                                                             // Fuzzy search technique
                FuzzySearch fuzzy;                                                            // Create FuzzySearch algorithm instance
                FuzzySearchResult fuzzy_result = fuzzy.search(sequence, pattern, max_errors); // Execute fuzzy search with edit distance
                result.positions = fuzzy_result.positions;                                   // Store all match positions found
                result.count = fuzzy_result.count;                                            // Store total number of matches
                result.success = true;                                                        // Mark execution as successful
                break;
            }
            case FUZZY_HAMMING: {                                                           // Fuzzy Hamming technique
                FuzzySearch fuzzy;                                                            // Create FuzzySearch algorithm instance
                FuzzySearchResult fuzzy_result = fuzzy.searchHamming(sequence, pattern, max_errors); // Execute Hamming distance search
                result.positions = fuzzy_result.positions;                                   // Store all match positions found
                result.count = fuzzy_result.count;                                            // Store total number of matches
                result.success = true;                                                        // Mark execution as successful
                break;
            }
            default:                                                                         // Unknown or unsupported technique
                result.success = false;                                                       // Mark execution as failed
                break;
        }
    } catch (...) {                                                                           // Catch any exceptions during execution
        result.success = false;                                                               // Mark as failed if exception occurs
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();                               // Record end time for this technique
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); // Calculate execution duration
    result.execution_time = duration.count();                                                 // Store execution time in microseconds
    
    return result;                                                                             // Return technique-specific result
}

void ConcurrentMultiSearch::combineResults(const std::vector<TechniqueResult>& results,
                                           MultiTechniqueResult& combined_result,
                                           const SearchConfig& config) {
    std::map<int, int> position_counts;                                                       // Map: position -> number of techniques that found it
    
    // Collect all positions from all techniques
    for (const auto& tech_result : results) {                                                 // Iterate through results from each search technique
        if (!tech_result.success) {                                                           // Check if technique execution succeeded
            continue;                                                                         // Skip failed techniques
        }
        
        combined_result.technique_positions[tech_result.technique_name] = tech_result.positions; // Store positions found by this technique
        combined_result.technique_counts[tech_result.technique_name] = tech_result.count;     // Store match count for this technique
        combined_result.technique_times[tech_result.technique_name] = tech_result.execution_time; // Store execution time for this technique
        
        // Count occurrences of each position across all techniques
        for (int pos : tech_result.positions) {                                               // Iterate through each position found by this technique
            position_counts[pos]++;                                                            // Increment count: how many techniques found this position
        }
    }
    
    // Find consensus positions (found by multiple techniques)
    if (config.require_consensus) {                                                           // Check if configuration requires consensus
        combined_result.consensus_positions = findConsensusPositions(                         // Find positions found by multiple techniques
                combined_result.technique_positions, config.min_consensus_count);             // Using minimum consensus count threshold
    } else {
        // Find positions found by at least min_consensus_count techniques
        for (const auto& pos_count : position_counts) {                                       // Iterate through position occurrence counts
            if (pos_count.second >= config.min_consensus_count) {                            // Check if position found by enough techniques
                combined_result.consensus_positions.push_back(pos_count.first);              // Add position to consensus list
            }
        }
    }
    
    // Find unique positions (found by only one technique)
    combined_result.unique_positions = findUniquePositions(                                  // Find positions found by exactly one technique
            combined_result.technique_positions);                                             // From all technique position maps
    
    // Calculate total unique matches
    std::set<int> all_positions;                                                              // Use set to automatically track unique positions
    for (const auto& tech_pos : combined_result.technique_positions) {                         // Iterate through all techniques' position lists
        for (int pos : tech_pos.second) {                                                     // Iterate through positions found by this technique
            all_positions.insert(pos);                                                         // Insert into set (duplicates automatically removed)
        }
    }
    combined_result.total_unique_matches = static_cast<int>(all_positions.size());            // Count total unique match positions
}

std::vector<int> ConcurrentMultiSearch::findConsensusPositions(
        const std::map<std::string, std::vector<int>>& technique_positions,
        int min_consensus) {
    std::map<int, int> position_counts;                                                       // Map: position -> number of techniques that found it
    
    // Count how many techniques found each position
    for (const auto& tech_pair : technique_positions) {                                      // Iterate through each technique's position list
        for (int pos : tech_pair.second) {                                                    // Iterate through positions found by this technique
            position_counts[pos]++;                                                           // Increment count: track how many techniques found this position
        }
    }
    
    std::vector<int> consensus;                                                               // Initialize vector for consensus positions
    for (const auto& pos_count : position_counts) {                                          // Iterate through position occurrence counts
        if (pos_count.second >= min_consensus) {                                            // Check if position found by enough techniques (meets threshold)
            consensus.push_back(pos_count.first);                                            // Add position to consensus list
        }
    }
    
    std::sort(consensus.begin(), consensus.end());                                            // Sort consensus positions in ascending order
    return consensus;                                                                         // Return sorted list of consensus positions
}

std::vector<int> ConcurrentMultiSearch::findUniquePositions(
        const std::map<std::string, std::vector<int>>& technique_positions) {
    std::map<int, int> position_counts;                                                       // Map: position -> number of techniques that found it
    
    // Count how many techniques found each position
    for (const auto& tech_pair : technique_positions) {                                      // Iterate through each technique's position list
        for (int pos : tech_pair.second) {                                                    // Iterate through positions found by this technique
            position_counts[pos]++;                                                           // Increment count: track how many techniques found this position
        }
    }
    
    std::vector<int> unique;                                                                  // Initialize vector for unique positions
    for (const auto& pos_count : position_counts) {                                          // Iterate through position occurrence counts
        if (pos_count.second == 1) {                                                         // Check if position found by exactly one technique
            unique.push_back(pos_count.first);                                                // Add position to unique list
        }
    }
    
    std::sort(unique.begin(), unique.end());                                                  // Sort unique positions in ascending order
    return unique;                                                                            // Return sorted list of unique positions
}

std::string ConcurrentMultiSearch::getTechniqueName(SearchTechnique technique) {
    switch (technique) {                                                                      // Switch on technique enumeration value
        case EXACT_MATCH:                                                                     // Exact match algorithm
            return "ExactMatch";                                                              // Return human-readable technique name
        case NAIVE_SEARCH:                                                                    // Naive brute-force search
            return "NaiveSearch";                                                            // Return human-readable technique name
        case RABIN_KARP:                                                                      // Rabin-Karp rolling hash
            return "RabinKarp";                                                              // Return human-readable technique name
        case KMP:                                                                             // Knuth-Morris-Pratt algorithm
            return "KMP";                                                                    // Return human-readable technique name
        case BOYER_MOORE:                                                                     // Boyer-Moore heuristic algorithm
            return "BoyerMoore";                                                             // Return human-readable technique name
        case FUZZY_SEARCH:                                                                    // Fuzzy search with edit distance
            return "FuzzySearch";                                                            // Return human-readable technique name
        case FUZZY_HAMMING:                                                                   // Fuzzy search with Hamming distance
            return "FuzzyHamming";                                                            // Return human-readable technique name
        case ALL_TECHNIQUES:                                                                  // All available techniques
            return "AllTechniques";                                                           // Return human-readable technique name
        default:                                                                              // Unknown or unsupported technique
            return "Unknown";                                                                 // Return default name for unknown techniques
    }
}

