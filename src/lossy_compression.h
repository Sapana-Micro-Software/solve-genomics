#ifndef LOSSY_COMPRESSION_H
#define LOSSY_COMPRESSION_H

#include "compression.h"
#include <string>
#include <map>
#include <vector>

/**
 * Lossy compression that approximates DNA sequences
 * Uses techniques like:
 * - Pattern approximation (replace similar patterns)
 * - Frequency-based compression (keep only frequent patterns)
 * - Truncation (remove low-information regions)
 */
class LossyCompression : public CompressionAlgorithm {
public:
    enum CompressionMode {
        FREQUENCY_BASED,    // Keep only frequent patterns
        PATTERN_APPROX,     // Approximate similar patterns
        TRUNCATION          // Remove low-information regions
    };
    
    LossyCompression(CompressionMode mode = FREQUENCY_BASED, 
                     double compression_target = 0.5);
    
    CompressionResult compress(const std::string& sequence) override;
    std::string decompress(const CompressionResult& compressed) override;
    std::string getName() const override { return "LossyCompression"; }
    bool isLossless() const override { return false; }
    
    /**
     * Set compression target ratio (0.0 to 1.0)
     */
    void setCompressionTarget(double target);
    
    /**
     * Frequency-based compression: keep only most frequent patterns
     */
    CompressionResult compressFrequencyBased(const std::string& sequence);
    
    /**
     * Pattern approximation: replace similar patterns with representatives
     */
    CompressionResult compressPatternApprox(const std::string& sequence);
    
    /**
     * Truncation: remove low-entropy regions
     */
    CompressionResult compressTruncation(const std::string& sequence);
    
private:
    CompressionMode mode_;
    double compression_target_;  // Target compression ratio (0.0 to 1.0)
    
    /**
     * Find similar patterns (within edit distance threshold)
     */
    std::map<std::string, std::string> findSimilarPatterns(
        const std::string& sequence,
        int pattern_length,
        int max_distance);
    
    /**
     * Calculate pattern frequencies
     */
    std::map<std::string, int> calculatePatternFrequencies(
        const std::string& sequence,
        int pattern_length);
    
    /**
     * Calculate local entropy for truncation
     */
    std::vector<double> calculateLocalEntropy(
        const std::string& sequence,
        int window_size);
};

#endif // LOSSY_COMPRESSION_H

