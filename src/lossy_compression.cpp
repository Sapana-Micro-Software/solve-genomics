#include "lossy_compression.h"
#include "sequence_generator.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <cstdlib>

LossyCompression::LossyCompression(CompressionMode mode, double compression_target)
    : mode_(mode), compression_target_(compression_target) {
    if (compression_target_ < 0.0) compression_target_ = 0.0;
    if (compression_target_ > 1.0) compression_target_ = 1.0;
}

CompressionResult LossyCompression::compress(const std::string& sequence) {
    switch (mode_) {
        case FREQUENCY_BASED:
            return compressFrequencyBased(sequence);
        case PATTERN_APPROX:
            return compressPatternApprox(sequence);
        case TRUNCATION:
            return compressTruncation(sequence);
        default:
            return compressFrequencyBased(sequence);
    }
}

std::string LossyCompression::decompress(const CompressionResult& compressed) {
    // For lossy compression, decompression returns approximation
    std::string decompressed = compressed.compressed_data;
    
    // Expand using dictionary/association list
    for (const auto& entry : compressed.dictionary) {
        const std::string& symbol = entry.first;
        const std::string& sequence = entry.second;
        
        size_t pos = 0;
        while ((pos = decompressed.find(symbol, pos)) != std::string::npos) {
            decompressed.replace(pos, symbol.length(), sequence);
            pos += sequence.length();
        }
    }
    
    return decompressed;
}

void LossyCompression::setCompressionTarget(double target) {
    compression_target_ = target;
    if (compression_target_ < 0.0) compression_target_ = 0.0;
    if (compression_target_ > 1.0) compression_target_ = 1.0;
}

CompressionResult LossyCompression::compressFrequencyBased(const std::string& sequence) {
    CompressionResult result;
    result.is_lossless = false;
    result.original_size = sequence.length();
    
    if (sequence.empty()) {
        result.compressed_data = "";
        result.compressed_size = 0;
        result.compression_ratio = 1.0;
        return result;
    }
    
    // Calculate pattern frequencies (using 4-base patterns)
    int pattern_length = 4;
    std::map<std::string, int> frequencies = 
        calculatePatternFrequencies(sequence, pattern_length);
    
    // Keep only top patterns that cover target compression ratio
    std::vector<std::pair<std::string, int>> sorted_patterns;
    for (const auto& freq : frequencies) {
        sorted_patterns.push_back(freq);
    }
    
    std::sort(sorted_patterns.begin(), sorted_patterns.end(),
        [](const std::pair<std::string, int>& a, 
           const std::pair<std::string, int>& b) {
            return a.second > b.second;
        });
    
    // Select patterns to keep
    size_t target_size = static_cast<size_t>(sequence.length() * compression_target_);
    std::map<std::string, std::string> pattern_to_symbol;
    std::string compressed = sequence;
    int symbol_num = 1;
    size_t current_size = sequence.length();
    
    for (const auto& pattern_entry : sorted_patterns) {
        if (current_size <= target_size) break;
        
        const std::string& pattern = pattern_entry.first;
        std::string symbol = "L" + std::to_string(symbol_num++);
        
        pattern_to_symbol[pattern] = symbol;
        result.dictionary[symbol] = pattern;
        
        // Replace pattern with symbol
        size_t pos = 0;
        while ((pos = compressed.find(pattern, pos)) != std::string::npos) {
            compressed.replace(pos, pattern.length(), symbol);
            current_size = current_size - pattern.length() + symbol.length();
            pos += symbol.length();
        }
    }
    
    result.compressed_data = compressed;
    result.compressed_size = compressed.length();
    result.compression_ratio = calculateCompressionRatio(
        result.original_size, result.compressed_size);
    
    return result;
}

CompressionResult LossyCompression::compressPatternApprox(const std::string& sequence) {
    CompressionResult result;
    result.is_lossless = false;
    result.original_size = sequence.length();
    
    if (sequence.empty()) {
        result.compressed_data = "";
        result.compressed_size = 0;
        result.compression_ratio = 1.0;
        return result;
    }
    
    // Find similar patterns and replace with representatives
    int pattern_length = 4;
    int max_distance = 1;  // Allow 1 edit distance for approximation
    
    std::map<std::string, std::string> similar_patterns = 
        findSimilarPatterns(sequence, pattern_length, max_distance);
    
    std::string compressed = sequence;
    int symbol_num = 1;
    
    // Replace similar patterns with representative
    for (const auto& entry : similar_patterns) {
        const std::string& pattern = entry.first;
        const std::string& representative = entry.second;
        
        if (pattern != representative) {
            std::string symbol = "A" + std::to_string(symbol_num++);
            result.dictionary[symbol] = representative;
            
            // Replace pattern with symbol (approximation)
            size_t pos = 0;
            while ((pos = compressed.find(pattern, pos)) != std::string::npos) {
                compressed.replace(pos, pattern.length(), symbol);
                pos += symbol.length();
            }
        }
    }
    
    result.compressed_data = compressed;
    result.compressed_size = compressed.length();
    result.compression_ratio = calculateCompressionRatio(
        result.original_size, result.compressed_size);
    
    return result;
}

CompressionResult LossyCompression::compressTruncation(const std::string& sequence) {
    CompressionResult result;
    result.is_lossless = false;
    result.original_size = sequence.length();
    
    if (sequence.empty()) {
        result.compressed_data = "";
        result.compressed_size = 0;
        result.compression_ratio = 1.0;
        return result;
    }
    
    // Calculate local entropy
    int window_size = 10;
    std::vector<double> local_entropy = calculateLocalEntropy(sequence, window_size);
    
    // Keep only high-entropy regions
    double entropy_threshold = 1.0;  // Threshold for keeping regions
    std::string compressed;
    
    for (size_t i = 0; i < sequence.length(); ++i) {
        if (local_entropy[i] >= entropy_threshold) {
            compressed += sequence[i];
        }
        // Low entropy regions are truncated (lossy)
    }
    
    result.compressed_data = compressed;
    result.compressed_size = compressed.length();
    result.compression_ratio = calculateCompressionRatio(
        result.original_size, result.compressed_size);
    
    return result;
}

std::map<std::string, std::string> LossyCompression::findSimilarPatterns(
    const std::string& sequence,
    int pattern_length,
    int max_distance) {
    
    std::map<std::string, std::string> similar_map;
    std::map<std::string, std::vector<std::string>> pattern_groups;
    
    // Find all patterns
    for (size_t i = 0; i <= sequence.length() - pattern_length; ++i) {
        std::string pattern = sequence.substr(i, pattern_length);
        
        // Find similar patterns (within edit distance)
        bool found_group = false;
        for (auto& group : pattern_groups) {
            const std::string& representative = group.first;
            
            // Calculate edit distance (simplified)
            int distance = 0;
            for (size_t j = 0; j < pattern.length() && j < representative.length(); ++j) {
                if (std::toupper(pattern[j]) != std::toupper(representative[j])) {
                    distance++;
                }
            }
            distance += std::abs(static_cast<int>(pattern.length() - representative.length()));
            
            if (distance <= max_distance) {
                group.second.push_back(pattern);
                similar_map[pattern] = representative;
                found_group = true;
                break;
            }
        }
        
        if (!found_group) {
            pattern_groups[pattern].push_back(pattern);
            similar_map[pattern] = pattern;  // Self-representative
        }
    }
    
    return similar_map;
}

std::map<std::string, int> LossyCompression::calculatePatternFrequencies(
    const std::string& sequence,
    int pattern_length) {
    
    std::map<std::string, int> frequencies;
    
    for (size_t i = 0; i <= sequence.length() - pattern_length; ++i) {
        std::string pattern = sequence.substr(i, pattern_length);
        frequencies[pattern]++;
    }
    
    return frequencies;
}

std::vector<double> LossyCompression::calculateLocalEntropy(
    const std::string& sequence,
    int window_size) {
    
    std::vector<double> local_entropy(sequence.length(), 0.0);
    
    for (size_t i = 0; i < sequence.length(); ++i) {
        size_t start = (i >= static_cast<size_t>(window_size)) ? 
                       (i - window_size) : 0;
        size_t end = std::min(i + window_size, sequence.length());
        
        std::string window = sequence.substr(start, end - start);
        local_entropy[i] = SequenceGenerator::calculateEntropy(window);
    }
    
    return local_entropy;
}

