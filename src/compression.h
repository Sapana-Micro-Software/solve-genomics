#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <string>
#include <vector>
#include <memory>
#include <map>

/**
 * Compression result structure
 */
struct CompressionResult {
    std::string compressed_data;      // Compressed representation
    std::vector<std::string> grammar; // Grammar rules (for grammar-based compression)
    std::map<std::string, std::string> dictionary; // Dictionary/association list
    size_t original_size;             // Original sequence size in bytes
    size_t compressed_size;           // Compressed size in bytes
    double compression_ratio;         // compressed_size / original_size
    bool is_lossless;                 // Whether compression is lossless
    
    CompressionResult() : original_size(0), compressed_size(0), 
                         compression_ratio(1.0), is_lossless(true) {}
};

/**
 * Base class for compression algorithms
 */
class CompressionAlgorithm {
public:
    virtual ~CompressionAlgorithm() = default;
    
    /**
     * Compress a DNA sequence
     * @param sequence Input DNA sequence
     * @return CompressionResult with compressed representation
     */
    virtual CompressionResult compress(const std::string& sequence) = 0;
    
    /**
     * Decompress compressed data back to original sequence
     * @param compressed Compressed data
     * @return Original sequence (or approximation if lossy)
     */
    virtual std::string decompress(const CompressionResult& compressed) = 0;
    
    /**
     * Get compression algorithm name
     */
    virtual std::string getName() const = 0;
    
    /**
     * Check if algorithm is lossless
     */
    virtual bool isLossless() const = 0;
    
    /**
     * Calculate compression ratio
     */
    static double calculateCompressionRatio(size_t original, size_t compressed) {
        if (original == 0) return 1.0;
        return static_cast<double>(compressed) / original;
    }
};

#endif // COMPRESSION_H

