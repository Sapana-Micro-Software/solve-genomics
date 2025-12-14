#ifndef PARALLEL_SEARCH_H
#define PARALLEL_SEARCH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <future>
#include <thread>
#include <mutex>
#include <atomic>

/**
 * Parallel, distributed, and concurrent search methods
 * Designed to scale to infinity for large patterns
 */
class ParallelSearch {
public:
    /**
     * Parallel search using multiple threads
     * @param sequence Sequence to search in
     * @param pattern Pattern to find
     * @param num_threads Number of threads
     * @return SearchResult with all matches
     */
    static SearchResult parallelSearch(const std::string& sequence,
                                      const std::string& pattern,
                                      int num_threads = std::thread::hardware_concurrency());
    
    /**
     * Distributed search simulation (chunk-based)
     * Divides sequence into chunks and processes in parallel
     * @param sequence Sequence to search
     * @param pattern Pattern to find
     * @param chunk_size Size of each chunk
     * @return SearchResult with all matches
     */
    static SearchResult distributedSearch(const std::string& sequence,
                                         const std::string& pattern,
                                         size_t chunk_size = 1000);
    
    /**
     * Concurrent search with async futures
     * @param sequence Sequence to search
     * @param pattern Pattern to find
     * @param num_tasks Number of concurrent tasks
     * @return SearchResult with all matches
     */
    static SearchResult concurrentSearch(const std::string& sequence,
                                        const std::string& pattern,
                                        int num_tasks = 4);
    
    /**
     * Map-Reduce style search
     * Maps pattern matching to chunks, reduces results
     * @param sequence Sequence to search
     * @param pattern Pattern to find
     * @param num_mappers Number of mapper threads
     * @return SearchResult with all matches
     */
    static SearchResult mapReduceSearch(const std::string& sequence,
                                       const std::string& pattern,
                                       int num_mappers = 4);
    
    /**
     * Pipeline search (producer-consumer pattern)
     * Producer generates chunks, consumers process them
     * @param sequence Sequence to search
     * @param pattern Pattern to find
     * @param num_consumers Number of consumer threads
     * @return SearchResult with all matches
     */
    static SearchResult pipelineSearch(const std::string& sequence,
                                      const std::string& pattern,
                                      int num_consumers = 4);
    
    /**
     * Work-stealing search (load balancing)
     * Threads steal work from others when idle
     * @param sequence Sequence to search
     * @param pattern Pattern to find
     * @param num_threads Number of threads
     * @return SearchResult with all matches
     */
    static SearchResult workStealingSearch(const std::string& sequence,
                                          const std::string& pattern,
                                          int num_threads = 4);
    
private:
    /**
     * Search in a chunk of sequence
     */
    static SearchResult searchChunk(const std::string& sequence,
                                   size_t start_pos,
                                   size_t end_pos,
                                   const std::string& pattern);
    
    /**
     * Merge search results from multiple chunks
     */
    static SearchResult mergeResults(const std::vector<SearchResult>& results);
};

#endif // PARALLEL_SEARCH_H

