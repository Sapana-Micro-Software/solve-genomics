#ifndef FUZZY_SEARCH_H
#define FUZZY_SEARCH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>

/**
 * Structure to hold fuzzy search results
 */
struct FuzzySearchResult {
    std::vector<int> positions;          // Positions where pattern was found
    std::vector<int> distances;          // Edit distance at each position
    int count;                           // Number of matches found
    int max_distance;                    // Maximum allowed distance used
    
    FuzzySearchResult() : count(0), max_distance(0) {}
};

/**
 * Fuzzy Search Algorithm
 * Finds approximate matches of a pattern in a sequence using edit distance
 * Supports both Hamming distance (for same-length strings) and 
 * Levenshtein distance (edit distance with insertions/deletions)
 */
class FuzzySearch {
public:
    /**
     * Search for pattern in sequence with maximum edit distance
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @param max_distance Maximum allowed edit distance (default: 0, exact match)
     * @return FuzzySearchResult containing all matches within distance threshold
     */
    FuzzySearchResult search(const std::string& sequence, 
                            const std::string& pattern, 
                            int max_distance = 0);
    
    /**
     * Search using Hamming distance (only substitutions, same length)
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @param max_distance Maximum allowed Hamming distance
     * @return FuzzySearchResult containing all matches
     */
    FuzzySearchResult searchHamming(const std::string& sequence,
                                   const std::string& pattern,
                                   int max_distance = 0);
    
    /**
     * Check if pattern exists within distance threshold
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @param max_distance Maximum allowed edit distance
     * @return true if pattern found within distance, false otherwise
     */
    bool contains(const std::string& sequence, 
                 const std::string& pattern, 
                 int max_distance = 0);
    
    /**
     * Find first occurrence within distance threshold
     * @param sequence The DNA sequence to search in
     * @param pattern The pattern to search for
     * @param max_distance Maximum allowed edit distance
     * @return Position of first match, or -1 if not found
     */
    int findFirst(const std::string& sequence,
                 const std::string& pattern,
                 int max_distance = 0);
    
    /**
     * Calculate edit distance (Levenshtein distance) between two strings
     * @param str1 First string
     * @param str2 Second string
     * @return Edit distance (minimum operations: insert, delete, substitute)
     */
    static int editDistance(const std::string& str1, const std::string& str2);
    
    /**
     * Calculate Hamming distance between two strings of same length
     * @param str1 First string
     * @param str2 Second string
     * @return Hamming distance (number of differing characters), or -1 if lengths differ
     */
    static int hammingDistance(const std::string& str1, const std::string& str2);
    
private:
    /**
     * Calculate edit distance with early termination if exceeds threshold
     * @param str1 First string
     * @param str2 Second string
     * @param max_distance Maximum distance to consider
     * @return Edit distance, or max_distance+1 if exceeds threshold
     */
    static int editDistanceBounded(const std::string& str1, 
                                   const std::string& str2,
                                   int max_distance);
};

#endif // FUZZY_SEARCH_H

