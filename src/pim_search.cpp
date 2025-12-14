#include "pim_search.h"
#include <algorithm>
#include <functional>
#include <cstring>

PIMSearch::PIMSearch() {
}

void PIMSearch::buildIndex(const std::vector<std::string>& sequences) {
    indexed_sequences_ = sequences;
    buildSuffixArray();
}

std::vector<size_t> PIMSearch::search(const std::string& pattern) {
    if (indexed_sequences_.empty()) {
        return std::vector<size_t>();
    }
    
    // Use first indexed sequence for search
    return vectorizedSearch(indexed_sequences_[0], pattern);
}

std::vector<std::vector<size_t>> PIMSearch::batchSearch(const std::vector<std::string>& patterns) {
    std::vector<std::vector<size_t>> results;
    
    for (const std::string& pattern : patterns) {
        results.push_back(search(pattern));
    }
    
    return results;
}

std::vector<size_t> PIMSearch::vectorizedSearch(const std::string& sequence, const std::string& pattern) {
    std::vector<size_t> positions;
    
    if (pattern.empty() || pattern.length() > sequence.length()) {
        return positions;
    }
    
    // Vectorized approach: process multiple positions at once
    // Simulate SIMD operations by processing in chunks
    const int chunk_size = 8;  // Process 8 characters at a time
    
    for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {
        bool match = true;
        
        // Process in chunks for better cache locality
        size_t j = 0;
        while (j < pattern.length() && match) {
            size_t chunk_end = std::min(j + chunk_size, pattern.length());
            
            // Compare chunk
            for (size_t k = j; k < chunk_end; ++k) {
                if (std::toupper(sequence[i + k]) != std::toupper(pattern[k])) {
                    match = false;
                    break;
                }
            }
            
            j = chunk_end;
        }
        
        if (match) {
            positions.push_back(i);
        }
    }
    
    return positions;
}

std::vector<size_t> PIMSearch::cacheOptimizedSearch(const std::string& pattern) {
    if (indexed_sequences_.empty()) {
        return std::vector<size_t>();
    }
    
    const std::string& sequence = indexed_sequences_[0];
    std::vector<size_t> positions;
    
    // Cache-optimized: process sequence in blocks that fit in cache
    const size_t cache_block_size = 64;  // Assume 64-byte cache line
    
    // Precompute pattern hash for quick rejection
    size_t pattern_hash = hashPattern(pattern);
    
    for (size_t block_start = 0; block_start < sequence.length(); block_start += cache_block_size) {
        size_t block_end = std::min(block_start + cache_block_size + pattern.length() - 1, 
                                   sequence.length());
        
        // Search within this cache block
        for (size_t i = block_start; i <= block_end - pattern.length(); ++i) {
            // Quick hash check before full comparison
            std::string window = sequence.substr(i, pattern.length());
            size_t window_hash = hashPattern(window);
            
            if (window_hash == pattern_hash) {
                // Hash match - verify with full comparison
                bool match = true;
                for (size_t j = 0; j < pattern.length(); ++j) {
                    if (std::toupper(sequence[i + j]) != std::toupper(pattern[j])) {
                        match = false;
                        break;
                    }
                }
                
                if (match) {
                    positions.push_back(i);
                }
            }
        }
    }
    
    return positions;
}

std::vector<size_t> PIMSearch::parallelSearch(const std::string& pattern, int num_threads) {
    if (indexed_sequences_.empty()) {
        return std::vector<size_t>();
    }
    
    const std::string& sequence = indexed_sequences_[0];
    std::vector<size_t> positions;
    
    // Simulate parallel processing by dividing sequence into chunks
    size_t chunk_size = sequence.length() / num_threads;
    
    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = (t == num_threads - 1) ? sequence.length() : (t + 1) * chunk_size;
        
        // Search in this chunk
        for (size_t i = start; i <= end - pattern.length(); ++i) {
            bool match = true;
            for (size_t j = 0; j < pattern.length(); ++j) {
                if (std::toupper(sequence[i + j]) != std::toupper(pattern[j])) {
                    match = false;
                    break;
                }
            }
            
            if (match) {
                positions.push_back(i);
            }
        }
    }
    
    // Sort positions (in real parallel implementation, would merge results)
    std::sort(positions.begin(), positions.end());
    
    return positions;
}

size_t PIMSearch::hashPattern(const std::string& pattern) {
    // Simple hash function for quick rejection
    std::hash<std::string> hasher;
    return hasher(pattern);
}

bool PIMSearch::vectorizedCompare(const std::string& str1, size_t pos1,
                                 const std::string& str2, size_t pos2,
                                 size_t length) {
    // Vectorized comparison (simulated)
    for (size_t i = 0; i < length; ++i) {
        if (std::toupper(str1[pos1 + i]) != std::toupper(str2[pos2 + i])) {
            return false;
        }
    }
    return true;
}

void PIMSearch::buildSuffixArray() {
    if (indexed_sequences_.empty()) {
        return;
    }
    
    const std::string& sequence = indexed_sequences_[0];
    suffixes_.clear();
    suffix_array_.clear();
    
    // Generate all suffixes
    for (size_t i = 0; i < sequence.length(); ++i) {
        suffixes_.push_back(sequence.substr(i));
        suffix_array_.push_back(i);
    }
    
    // Sort suffixes (in real implementation, would use efficient suffix array construction)
    std::sort(suffix_array_.begin(), suffix_array_.end(),
        [&sequence](size_t a, size_t b) {
            return sequence.substr(a) < sequence.substr(b);
        });
}

