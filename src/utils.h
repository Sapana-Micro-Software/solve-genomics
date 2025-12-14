#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

/**
 * Utility functions for DNA sequence processing
 */
namespace DNAUtils {

    /**
     * Validate that a sequence contains only valid DNA bases (A, T, G, C)
     * @param sequence DNA sequence to validate
     * @return true if valid, false otherwise
     */
    bool isValidDNA(const std::string& sequence);
    
    /**
     * Convert sequence to uppercase
     * @param sequence Input sequence
     * @return Uppercase sequence
     */
    std::string toUpperCase(const std::string& sequence);
    
    /**
     * Remove whitespace from sequence
     * @param sequence Input sequence
     * @return Sequence without whitespace
     */
    std::string removeWhitespace(const std::string& sequence);
    
    /**
     * Clean and validate DNA sequence
     * @param sequence Input sequence
     * @return Cleaned sequence, or empty string if invalid
     */
    std::string cleanSequence(const std::string& sequence);
    
    /**
     * Read DNA sequence from file
     * @param filename Path to file
     * @return Sequence read from file, or empty string if error
     */
    std::string readSequenceFromFile(const std::string& filename);
    
    /**
     * Format alignment for display
     * @param seq1 First aligned sequence
     * @param seq2 Second aligned sequence
     * @param line_length Maximum characters per line
     * @return Formatted string showing alignment
     */
    std::string formatAlignment(const std::string& seq1, const std::string& seq2, int line_length = 80);
    
    /**
     * Calculate percentage identity between two aligned sequences
     * @param aligned_seq1 First aligned sequence
     * @param aligned_seq2 Second aligned sequence
     * @return Percentage identity (0-100)
     */
    double calculateIdentity(const std::string& aligned_seq1, const std::string& aligned_seq2);
}

#endif // UTILS_H

