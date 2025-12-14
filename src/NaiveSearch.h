#ifndef NAIVE_SEARCH_H
#define NAIVE_SEARCH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>

/**
 * Naive String Matching Algorithm
 * Brute-force pattern matching that checks all possible positions
 */
class NaiveSearch {
public:
    /**
     * Search for pattern in sequence using naive algorithm
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @return SearchResult containing all positions where pattern was found
     */
    SearchResult search(const std::string& sequence, const std::string& pattern);
    
    /**
     * Find first occurrence of pattern
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @return Position of first match, or -1 if not found
     */
    int findFirst(const std::string& sequence, const std::string& pattern);
};

#endif // NAIVE_SEARCH_H

