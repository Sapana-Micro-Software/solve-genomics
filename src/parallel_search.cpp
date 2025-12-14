#include "parallel_search.h"
#include "ExactMatch.h"
#include <algorithm>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <future>

SearchResult ParallelSearch::parallelSearch(const std::string& sequence,
                                           const std::string& pattern,
                                           int num_threads) {
    if (sequence.empty() || pattern.empty() || pattern.length() > sequence.length()) {
        return SearchResult();
    }
    
    if (num_threads <= 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 4;
    }
    
    size_t seq_len = sequence.length();
    size_t chunk_size = (seq_len + num_threads - 1) / num_threads;
    
    std::vector<std::future<SearchResult>> futures;
    
    // Launch parallel searches
    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size + pattern.length() - 1, seq_len);
        
        futures.push_back(std::async(std::launch::async, [&sequence, start, end, &pattern]() {
            return searchChunk(sequence, start, end, pattern);
        }));
    }
    
    // Collect results
    std::vector<SearchResult> results;
    for (auto& future : futures) {
        results.push_back(future.get());
    }
    
    return mergeResults(results);
}

SearchResult ParallelSearch::distributedSearch(const std::string& sequence,
                                              const std::string& pattern,
                                              size_t chunk_size) {
    if (sequence.empty() || pattern.empty()) {
        return SearchResult();
    }
    
    if (chunk_size == 0) {
        chunk_size = 1000;
    }
    
    std::vector<std::future<SearchResult>> futures;
    size_t seq_len = sequence.length();
    
    // Divide into chunks
    for (size_t start = 0; start < seq_len; start += chunk_size) {
        size_t end = std::min(start + chunk_size + pattern.length() - 1, seq_len);
        
        futures.push_back(std::async(std::launch::async, [&sequence, start, end, &pattern]() {
            return searchChunk(sequence, start, end, pattern);
        }));
    }
    
    // Collect and merge
    std::vector<SearchResult> results;
    for (auto& future : futures) {
        results.push_back(future.get());
    }
    
    return mergeResults(results);
}

SearchResult ParallelSearch::concurrentSearch(const std::string& sequence,
                                             const std::string& pattern,
                                             int num_tasks) {
    if (num_tasks <= 0) {
        num_tasks = 4;
    }
    
    return distributedSearch(sequence, pattern, sequence.length() / num_tasks);
}

SearchResult ParallelSearch::mapReduceSearch(const std::string& sequence,
                                            const std::string& pattern,
                                            int num_mappers) {
    if (num_mappers <= 0) {
        num_mappers = 4;
    }
    
    size_t seq_len = sequence.length();
    size_t chunk_size = (seq_len + num_mappers - 1) / num_mappers;
    
    // Map phase: parallel search in chunks
    std::vector<std::future<SearchResult>> map_futures;
    
    for (int i = 0; i < num_mappers; ++i) {
        size_t start = i * chunk_size;
        size_t end = std::min(start + chunk_size + pattern.length() - 1, seq_len);
        
        map_futures.push_back(std::async(std::launch::async, [&sequence, start, end, &pattern]() {
            return searchChunk(sequence, start, end, pattern);
        }));
    }
    
    // Reduce phase: merge all results
    std::vector<SearchResult> map_results;
    for (auto& future : map_futures) {
        map_results.push_back(future.get());
    }
    
    return mergeResults(map_results);
}

SearchResult ParallelSearch::pipelineSearch(const std::string& sequence,
                                           const std::string& pattern,
                                           int num_consumers) {
    if (num_consumers <= 0) {
        num_consumers = 4;
    }
    
    std::queue<size_t> chunk_queue;
    std::mutex queue_mutex;
    std::condition_variable queue_cv;
    std::vector<SearchResult> results;
    std::mutex results_mutex;
    bool done = false;
    
    size_t chunk_size = 1000;
    size_t seq_len = sequence.length();
    size_t current_chunk = 0;
    
    // Consumer threads
    std::vector<std::thread> consumers;
    
    for (int i = 0; i < num_consumers; ++i) {
        consumers.push_back(std::thread([&]() {
            while (true) {
                size_t chunk_start;
                
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    queue_cv.wait(lock, [&]() { return !chunk_queue.empty() || done; });
                    
                    if (chunk_queue.empty() && done) {
                        break;
                    }
                    
                    if (!chunk_queue.empty()) {
                        chunk_start = chunk_queue.front();
                        chunk_queue.pop();
                    } else {
                        continue;
                    }
                }
                
                // Process chunk
                size_t chunk_end = std::min(chunk_start + chunk_size + pattern.length() - 1, seq_len);
                SearchResult chunk_result = searchChunk(sequence, chunk_start, chunk_end, pattern);
                
                // Add to results
                {
                    std::lock_guard<std::mutex> lock(results_mutex);
                    results.push_back(chunk_result);
                }
            }
        }));
    }
    
    // Producer: add chunks to queue
    while (current_chunk < seq_len) {
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            chunk_queue.push(current_chunk);
        }
        queue_cv.notify_one();
        current_chunk += chunk_size;
    }
    
    // Signal done
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        done = true;
    }
    queue_cv.notify_all();
    
    // Wait for consumers
    for (auto& consumer : consumers) {
        consumer.join();
    }
    
    return mergeResults(results);
}

SearchResult ParallelSearch::workStealingSearch(const std::string& sequence,
                                               const std::string& pattern,
                                               int num_threads) {
    if (num_threads <= 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 4;
    }
    
    size_t seq_len = sequence.length();
    size_t total_chunks = (seq_len + 999) / 1000;  // ~1000 base chunks
    
    std::atomic<size_t> next_chunk(0);
    std::vector<SearchResult> thread_results(num_threads);
    std::vector<std::thread> threads;
    
    // Worker threads with work stealing
    for (int t = 0; t < num_threads; ++t) {
        threads.push_back(std::thread([&, t]() {
            while (true) {
                size_t chunk_idx = next_chunk.fetch_add(1);
                if (chunk_idx >= total_chunks) {
                    break;
                }
                
                size_t chunk_start = chunk_idx * 1000;
                size_t chunk_end = std::min(chunk_start + 1000 + pattern.length() - 1, seq_len);
                
                SearchResult chunk_result = searchChunk(sequence, chunk_start, chunk_end, pattern);
                
                // Merge into thread's result
                std::vector<SearchResult> to_merge = {thread_results[t], chunk_result};
                thread_results[t] = mergeResults(to_merge);
            }
        }));
    }
    
    // Wait for all threads
    for (auto& thread : threads) {
        thread.join();
    }
    
    return mergeResults(thread_results);
}

SearchResult ParallelSearch::searchChunk(const std::string& sequence,
                                        size_t start_pos,
                                        size_t end_pos,
                                        const std::string& pattern) {
    if (end_pos > sequence.length()) {
        end_pos = sequence.length();
    }
    
    if (start_pos >= end_pos || pattern.empty()) {
        return SearchResult();
    }
    
    std::string chunk = sequence.substr(start_pos, end_pos - start_pos);
    
    ExactMatch matcher;
    SearchResult result = matcher.search(chunk, pattern);
    
    // Adjust positions to global coordinates
    for (int& pos : result.positions) {
        pos += static_cast<int>(start_pos);
    }
    
    return result;
}

SearchResult ParallelSearch::mergeResults(const std::vector<SearchResult>& results) {
    SearchResult merged;
    
    for (const auto& result : results) {
        merged.count += result.count;
        merged.positions.insert(merged.positions.end(),
                               result.positions.begin(),
                               result.positions.end());
    }
    
    // Sort and remove duplicates
    std::sort(merged.positions.begin(), merged.positions.end());
    auto last = std::unique(merged.positions.begin(), merged.positions.end());
    merged.positions.erase(last, merged.positions.end());
    merged.count = merged.positions.size();
    
    return merged;
}

