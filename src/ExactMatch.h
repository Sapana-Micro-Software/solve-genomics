#ifndef EXACT_MATCH_H
#define EXACT_MATCH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>

/**
 * Exact Match Search Algorithm
 * Finds all exact occurrences of a pattern in a sequence
 */
class ExactMatch {
public:
    /**
     * Search for exact occurrences of pattern in sequence
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @return SearchResult containing all positions where pattern was found
     */
    SearchResult search(const std::string& sequence, const std::string& pattern);
    
    /**
     * Check if pattern exists in sequence
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @return true if pattern found, false otherwise
     */
    bool contains(const std::string& sequence, const std::string& pattern);
    
    /**
     * Find first occurrence of pattern in sequence
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @return Position of first match, or -1 if not found
     */
    int findFirst(const std::string& sequence, const std::string& pattern);
};

#endif // EXACT_MATCH_H

