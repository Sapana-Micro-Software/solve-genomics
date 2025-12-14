#ifndef EDIT_DISTANCE_ALGORITHMS_H
#define EDIT_DISTANCE_ALGORITHMS_H

#include <string>
#include <vector>

/**
 * Collection of edit distance algorithms with different optimizations
 */
class EditDistanceAlgorithms {
public:
    /**
     * Standard Levenshtein distance (full DP matrix)
     * Time: O(n*m), Space: O(n*m)
     */
    static int levenshteinDistance(const std::string& str1, const std::string& str2);
    
    /**
     * Space-optimized Levenshtein (only two rows)
     * Time: O(n*m), Space: O(min(n,m))
     */
    static int levenshteinDistanceOptimized(const std::string& str1, const std::string& str2);
    
    /**
     * Bounded Levenshtein with early termination
     * Time: O(n*m*d), Space: O(min(n,m)) where d is max_distance
     */
    static int levenshteinDistanceBounded(const std::string& str1, 
                                         const std::string& str2,
                                         int max_distance);
    
    /**
     * Damerau-Levenshtein distance (includes transpositions)
     * Time: O(n*m), Space: O(n*m)
     */
    static int damerauLevenshteinDistance(const std::string& str1, const std::string& str2);
    
    /**
     * Hamming distance (only substitutions, same length)
     * Time: O(n), Space: O(1)
     */
    static int hammingDistance(const std::string& str1, const std::string& str2);
    
    /**
     * Jaro-Winkler distance (similarity measure, 0.0 to 1.0)
     * Time: O(n*m), Space: O(n*m)
     */
    static double jaroWinklerDistance(const std::string& str1, const std::string& str2);
    
    /**
     * Longest Common Subsequence (LCS) length
     * Time: O(n*m), Space: O(n*m)
     */
    static int longestCommonSubsequence(const std::string& str1, const std::string& str2);
    
    /**
     * Longest Common Substring length
     * Time: O(n*m), Space: O(n*m)
     */
    static int longestCommonSubstring(const std::string& str1, const std::string& str2);
    
    /**
     * Weighted edit distance with custom costs
     * @param str1 First string
     * @param str2 Second string
     * @param insert_cost Cost of insertion
     * @param delete_cost Cost of deletion
     * @param substitute_cost Cost of substitution
     * @return Weighted edit distance
     */
    static int weightedEditDistance(const std::string& str1,
                                   const std::string& str2,
                                   int insert_cost = 1,
                                   int delete_cost = 1,
                                   int substitute_cost = 1);
    
    /**
     * Calculate edit distance with DNA-specific scoring
     * Different costs for transitions (A<->G, T<->C) vs transversions
     */
    static int dnaEditDistance(const std::string& seq1, const std::string& seq2);
    
private:
    /**
     * Check if two DNA bases can transition (A<->G or T<->C)
     */
    static bool isTransition(char base1, char base2);
};

#endif // EDIT_DISTANCE_ALGORITHMS_H

