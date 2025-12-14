#ifndef ASSOCIATION_LIST_H
#define ASSOCIATION_LIST_H

#include "compression.h"
#include <map>
#include <string>
#include <vector>

/**
 * Association list for pattern matching on compressed sequences
 * Maps compressed symbols/patterns to their original sequences
 */
class AssociationList {
public:
    /**
     * Add association (symbol -> sequence)
     */
    void addAssociation(const std::string& symbol, const std::string& sequence);
    
    /**
     * Get sequence for symbol
     */
    std::string getSequence(const std::string& symbol) const;
    
    /**
     * Check if symbol exists
     */
    bool hasSymbol(const std::string& symbol) const;
    
    /**
     * Expand compressed sequence using association list
     */
    std::string expand(const std::string& compressed_sequence) const;
    
    /**
     * Build association list from grammar rules
     */
    void buildFromGrammar(const std::vector<std::string>& grammar);
    
    /**
     * Build association list from dictionary
     */
    void buildFromDictionary(const std::map<std::string, std::string>& dictionary);
    
    /**
     * Get all symbols
     */
    std::vector<std::string> getSymbols() const;
    
    /**
     * Clear all associations
     */
    void clear();
    
    /**
     * Get size
     */
    size_t size() const;
    
private:
    std::map<std::string, std::string> associations_;  // symbol -> sequence mapping
};

/**
 * Pattern matcher that works on compressed sequences using association lists
 */
class CompressedPatternMatcher {
public:
    CompressedPatternMatcher(const AssociationList& association_list);
    
    /**
     * Search for pattern in compressed sequence
     * @param compressed_sequence Compressed sequence
     * @param pattern Pattern to search for
     * @return Positions where pattern matches (in original sequence coordinates)
     */
    std::vector<size_t> search(const std::string& compressed_sequence, 
                              const std::string& pattern);
    
    /**
     * Search with fuzzy matching (edit distance)
     * @param compressed_sequence Compressed sequence
     * @param pattern Pattern to search for
     * @param max_distance Maximum edit distance
     * @return Positions and distances
     */
    std::vector<std::pair<size_t, int>> searchFuzzy(
        const std::string& compressed_sequence,
        const std::string& pattern,
        int max_distance);
    
private:
    const AssociationList& association_list_;
    
    /**
     * Expand compressed sequence to original
     */
    std::string expandSequence(const std::string& compressed) const;
    
    /**
     * Convert compressed position to original position
     */
    size_t convertPosition(const std::string& compressed, size_t compressed_pos) const;
    
    /**
     * Calculate edit distance (helper function)
     */
    int calculateEditDistance(const std::string& str1, 
                             const std::string& str2,
                             int max_distance) const;
};

#endif // ASSOCIATION_LIST_H

