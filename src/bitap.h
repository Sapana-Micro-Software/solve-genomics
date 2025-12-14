#ifndef BITAP_H
#define BITAP_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>

/**
 * Bitap Algorithm (Shift-Or) for Approximate String Matching
 * Uses bit-parallelism for efficient fuzzy matching
 * Time: O(n*m/w) where w is word size (typically 32 or 64)
 * Space: O(Σ) where Σ is alphabet size
 */
class Bitap {
public:
    Bitap(int max_errors = 2);
    
    /**
     * Search for pattern with allowed errors
     * @param text Text to search in
     * @param pattern Pattern to find
     * @param max_errors Maximum number of errors allowed
     * @return SearchResult with matches
     */
    SearchResult search(const std::string& text, 
                       const std::string& pattern,
                       int max_errors = -1);
    
    /**
     * Check if pattern exists with errors
     * @param text Text to search
     * @param pattern Pattern to find
     * @param max_errors Maximum errors
     * @return True if found
     */
    bool contains(const std::string& text, 
                 const std::string& pattern,
                 int max_errors = -1);
    
    /**
     * Find first occurrence
     * @param text Text to search
     * @param pattern Pattern to find
     * @param max_errors Maximum errors
     * @return Position or -1 if not found
     */
    int findFirst(const std::string& text,
                 const std::string& pattern,
                 int max_errors = -1);
    
private:
    int max_errors_;
    
    /**
     * Build character mask table
     * @param pattern Pattern to build mask for
     * @return Map from character to bitmask
     */
    std::map<char, unsigned long long> buildCharMask(const std::string& pattern);
    
    /**
     * Shift-Or algorithm implementation
     */
    std::vector<size_t> shiftOrSearch(const std::string& text,
                                     const std::string& pattern,
                                     int max_errors);
    
    /**
     * Calculate edit distance using bit-parallelism
     */
    int bitapEditDistance(const std::string& text, const std::string& pattern);
};

#endif // BITAP_H

