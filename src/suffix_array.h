#ifndef SUFFIX_ARRAY_H
#define SUFFIX_ARRAY_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <algorithm>

/**
 * Suffix Array for efficient substring search
 * Time: O(n log n) for construction, O(m log n) for search
 * Space: O(n)
 */
class SuffixArray {
public:
    SuffixArray();
    
    /**
     * Build suffix array from text
     * @param text Text to build array from
     */
    void build(const std::string& text);
    
    /**
     * Search for pattern using binary search
     * @param pattern Pattern to search for
     * @return SearchResult with all occurrences
     */
    SearchResult search(const std::string& pattern);
    
    /**
     * Find first occurrence
     * @param pattern Pattern to find
     * @return First position or -1 if not found
     */
    int findFirst(const std::string& pattern);
    
    /**
     * Find longest common prefix (LCP) array
     * @return LCP array
     */
    std::vector<int> buildLCP();
    
    /**
     * Find longest repeated substring
     * @return Longest repeated substring
     */
    std::string longestRepeatedSubstring();
    
private:
    std::string text_;
    std::vector<int> suffix_array_;  // Sorted suffix indices
    std::vector<int> rank_;           // Rank of each suffix
    
    /**
     * Compare two suffixes
     */
    int compareSuffix(int i, int j, int k);
    
    /**
     * Build suffix array using doubling method
     */
    void buildDoubling();
    
    /**
     * Binary search for pattern
     */
    std::pair<int, int> binarySearch(const std::string& pattern);
    
    /**
     * Compare pattern with suffix
     */
    int comparePatternSuffix(const std::string& pattern, int suffix_idx);
};

#endif // SUFFIX_ARRAY_H

