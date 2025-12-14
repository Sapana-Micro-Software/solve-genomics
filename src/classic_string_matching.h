#ifndef CLASSIC_STRING_MATCHING_H
#define CLASSIC_STRING_MATCHING_H

#include "SequenceAligner.h"
#include <string>
#include <vector>

/**
 * Classic string matching algorithms for DNA sequence search
 */
class ClassicStringMatching {
public:
    /**
     * Rabin-Karp Algorithm
     * Uses rolling hash for efficient pattern matching
     * Time: O(n+m) average, O(n*m) worst case
     * Space: O(1)
     */
    class RabinKarp {
    public:
        /**
         * Search for pattern in text using Rabin-Karp
         * @param text Text to search in
         * @param pattern Pattern to search for
         * @return Vector of positions where pattern found
         */
        static SearchResult search(const std::string& text, const std::string& pattern);
        
        /**
         * Search with custom hash base and modulus
         * @param text Text to search in
         * @param pattern Pattern to search for
         * @param base Hash base (default: 256)
         * @param modulus Hash modulus (default: 101)
         * @return Vector of positions
         */
        static SearchResult searchWithHash(const std::string& text,
                                          const std::string& pattern,
                                          int base = 256,
                                          int modulus = 101);
        
    private:
        /**
         * Calculate hash value for string
         */
        static int calculateHash(const std::string& str, int base, int modulus);
        
        /**
         * Calculate rolling hash
         */
        static int rollingHash(int old_hash, char old_char, char new_char,
                              int base, int modulus, int pattern_len);
    };
    
    /**
     * Knuth-Morris-Pratt (KMP) Algorithm
     * Uses failure function for efficient pattern matching
     * Time: O(n+m)
     * Space: O(m)
     */
    class KMP {
    public:
        /**
         * Search for pattern in text using KMP
         * @param text Text to search in
         * @param pattern Pattern to search for
         * @return Vector of positions where pattern found
         */
        static SearchResult search(const std::string& text, const std::string& pattern);
        
        /**
         * Find first occurrence
         * @param text Text to search in
         * @param pattern Pattern to search for
         * @return Position of first match, or -1 if not found
         */
        static int findFirst(const std::string& text, const std::string& pattern);
        
    private:
        /**
         * Build failure function (LPS array - Longest Proper Prefix which is also Suffix)
         * @param pattern Pattern to build failure function for
         * @return Failure function array
         */
        static std::vector<int> buildFailureFunction(const std::string& pattern);
    };
    
    /**
     * Boyer-Moore Algorithm
     * Uses bad character and good suffix heuristics
     * Time: O(n*m) worst case, O(n/m) best case
     * Space: O(m + alphabet_size)
     */
    class BoyerMoore {
    public:
        /**
         * Search for pattern in text using Boyer-Moore
         * @param text Text to search in
         * @param pattern Pattern to search for
         * @return Vector of positions where pattern found
         */
        static SearchResult search(const std::string& text, const std::string& pattern);
        
        /**
         * Find first occurrence
         * @param text Text to search in
         * @param pattern Pattern to search for
         * @return Position of first match, or -1 if not found
         */
        static int findFirst(const std::string& text, const std::string& pattern);
        
    private:
        /**
         * Build bad character table
         * Maps each character to its rightmost position in pattern
         * @param pattern Pattern to build table for
         * @return Bad character table
         */
        static std::vector<int> buildBadCharacterTable(const std::string& pattern);
        
        /**
         * Build good suffix table
         * @param pattern Pattern to build table for
         * @return Good suffix table
         */
        static std::vector<int> buildGoodSuffixTable(const std::string& pattern);
        
        /**
         * Build suffix array for good suffix heuristic
         */
        static std::vector<int> buildSuffixArray(const std::string& pattern);
    };
    
    /**
     * Compare all three algorithms on same input
     * @param text Text to search in
     * @param pattern Pattern to search for
     * @return Results from all three algorithms
     */
    struct ComparisonResult {
        SearchResult rabin_karp;
        SearchResult kmp;
        SearchResult boyer_moore;
    };
    
    static ComparisonResult compareAll(const std::string& text, const std::string& pattern);
};

#endif // CLASSIC_STRING_MATCHING_H

