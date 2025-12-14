#ifndef WU_MANBER_H
#define WU_MANBER_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

/**
 * Wu-Manber Algorithm for Multi-Pattern Matching
 * Uses shift tables and hash tables for efficient multi-pattern search
 * Time: O(n*m/B + r) where n is text length, m is pattern length, B is block size, r is matches
 * Space: O(Σ * m) where Σ is alphabet size
 */
class WuManber {
public:
    struct MultiPatternResult {
        std::map<std::string, std::vector<size_t>> matches;  // Pattern -> positions
        int total_matches;
        
        MultiPatternResult() : total_matches(0) {}
    };
    
    WuManber(int block_size = 2);
    
    /**
     * Build shift and hash tables from patterns
     * @param patterns Vector of patterns to search for
     */
    void buildTables(const std::vector<std::string>& patterns);
    
    /**
     * Search for all patterns in text
     * @param text Text to search in
     * @return MultiPatternResult with all matches
     */
    MultiPatternResult search(const std::string& text);
    
    /**
     * Search for single pattern (for compatibility)
     * @param text Text to search in
     * @param pattern Pattern to find
     * @return SearchResult with matches
     */
    SearchResult searchSingle(const std::string& text, const std::string& pattern);
    
private:
    int block_size_;
    std::vector<std::string> patterns_;
    
    // Shift table: maps block to minimum shift
    std::unordered_map<std::string, int> shift_table_;
    
    // Hash table: maps block to pattern indices
    std::unordered_map<std::string, std::vector<int>> hash_table_;
    
    // Prefix table for verification
    std::vector<std::string> prefix_table_;
    
    /**
     * Get block from position
     */
    std::string getBlock(const std::string& text, size_t pos);
    
    /**
     * Calculate shift value
     */
    int calculateShift(const std::string& block, int pattern_index);
    
    /**
     * Verify match
     */
    bool verifyMatch(const std::string& text, size_t pos, int pattern_index);
};

#endif // WU_MANBER_H

