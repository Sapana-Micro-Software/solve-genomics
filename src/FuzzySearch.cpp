#include "FuzzySearch.h"
#include "utils.h"
#include <algorithm>
#include <climits>
#include <cstdlib>

FuzzySearchResult FuzzySearch::search(const std::string& sequence, 
                                      const std::string& pattern, 
                                      int max_distance) {
    FuzzySearchResult result;                                                                  // Initialize result structure
    result.max_distance = max_distance;                                                       // Store maximum allowed distance
    
    if (pattern.empty() || sequence.empty() || pattern.length() > sequence.length() + max_distance) { // Validate input parameters
        return result;                                                                        // Return empty result if invalid
    }
    
    size_t seq_len = sequence.length();                                                       // Get sequence length
    size_t pat_len = pattern.length();                                                        // Get pattern length
    
    // Slide pattern over sequence
    for (size_t i = 0; i <= seq_len - pat_len + max_distance; ++i) {                          // Slide window allowing for edit distance
        // Extract substring (may extend beyond sequence for edit distance)
        size_t end_pos = std::min(i + pat_len + max_distance, seq_len);                       // Calculate end position (don't exceed sequence)
        std::string substring = sequence.substr(i, end_pos - i);                              // Extract substring to compare
        
        // Calculate edit distance
        int distance = editDistanceBounded(pattern, substring, max_distance);                 // Compute edit distance with bound
        
        if (distance <= max_distance) {                                                        // Check if within threshold
            result.positions.push_back(static_cast<int>(i));                                 // Record position of match
            result.distances.push_back(distance);                                             // Record edit distance
            result.count++;                                                                   // Increment match counter
        }
    }
    
    return result;                                                                            // Return all matches within distance threshold
}

FuzzySearchResult FuzzySearch::searchHamming(const std::string& sequence,
                                             const std::string& pattern,
                                             int max_distance) {
    FuzzySearchResult result;                                                                  // Initialize result structure
    result.max_distance = max_distance;                                                       // Store maximum allowed distance
    
    if (pattern.empty() || sequence.empty() || pattern.length() > sequence.length()) {        // Validate input parameters
        return result;                                                                        // Return empty result if invalid
    }
    
    size_t seq_len = sequence.length();                                                       // Get sequence length
    size_t pat_len = pattern.length();                                                        // Get pattern length
    
    // Hamming distance requires same length, so slide exact-length window
    for (size_t i = 0; i <= seq_len - pat_len; ++i) {                                         // Slide exact-length window through sequence
        std::string substring = sequence.substr(i, pat_len);                                  // Extract substring of pattern length
        int distance = hammingDistance(pattern, substring);                                  // Calculate Hamming distance (substitutions only)
        
        if (distance != -1 && distance <= max_distance) {                                    // Check if valid and within threshold
            result.positions.push_back(static_cast<int>(i));                                 // Record position of match
            result.distances.push_back(distance);                                             // Record Hamming distance
            result.count++;                                                                   // Increment match counter
        }
    }
    
    return result;                                                                            // Return all matches within Hamming distance threshold
}

bool FuzzySearch::contains(const std::string& sequence, 
                           const std::string& pattern, 
                           int max_distance) {
    return findFirst(sequence, pattern, max_distance) != -1;
}

int FuzzySearch::findFirst(const std::string& sequence,
                          const std::string& pattern,
                          int max_distance) {
    if (pattern.empty() || sequence.empty() || pattern.length() > sequence.length() + max_distance) {
        return -1;
    }
    
    size_t seq_len = sequence.length();
    size_t pat_len = pattern.length();
    
    for (size_t i = 0; i <= seq_len - pat_len + max_distance; ++i) {
        size_t end_pos = std::min(i + pat_len + max_distance, seq_len);
        std::string substring = sequence.substr(i, end_pos - i);
        
        int distance = editDistanceBounded(pattern, substring, max_distance);
        
        if (distance <= max_distance) {
            return static_cast<int>(i);
        }
    }
    
    return -1;
}

int FuzzySearch::editDistance(const std::string& str1, const std::string& str2) {
    int m = str1.length();
    int n = str2.length();
    
    // Create DP table
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    // Initialize first row and column
    for (int i = 0; i <= m; ++i) {
        dp[i][0] = i;
    }
    for (int j = 0; j <= n; ++j) {
        dp[0][j] = j;
    }
    
    // Fill the table
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            // Case-insensitive comparison
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                dp[i][j] = dp[i - 1][j - 1];  // Match
            } else {
                // Take minimum of three operations
                dp[i][j] = 1 + std::min({
                    dp[i - 1][j],      // Delete from str1
                    dp[i][j - 1],      // Insert into str1
                    dp[i - 1][j - 1]   // Substitute
                });
            }
        }
    }
    
    return dp[m][n];
}

int FuzzySearch::hammingDistance(const std::string& str1, const std::string& str2) {
    if (str1.length() != str2.length()) {
        return -1;  // Hamming distance only defined for same-length strings
    }
    
    int distance = 0;
    for (size_t i = 0; i < str1.length(); ++i) {
        if (std::toupper(str1[i]) != std::toupper(str2[i])) {
            distance++;
        }
    }
    
    return distance;
}

int FuzzySearch::editDistanceBounded(const std::string& str1, 
                                     const std::string& str2,
                                     int max_distance) {
    int m = str1.length();
    int n = str2.length();
    
    // Early termination: if length difference exceeds max_distance, 
    // edit distance must be at least that
    int len_diff = (m > n) ? (m - n) : (n - m);
    if (len_diff > max_distance) {
        return max_distance + 1;
    }
    
    // Use space-optimized DP (only need two rows)
    std::vector<int> prev(n + 1);
    std::vector<int> curr(n + 1);
    
    // Initialize first row
    for (int j = 0; j <= n; ++j) {
        prev[j] = j;
    }
    
    // Fill the table
    for (int i = 1; i <= m; ++i) {
        curr[0] = i;
        
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                curr[j] = prev[j - 1];  // Match
            } else {
                curr[j] = 1 + std::min({
                    prev[j],      // Delete
                    curr[j - 1],  // Insert
                    prev[j - 1]   // Substitute
                });
            }
        }
        
        // Early termination: if minimum in current row exceeds threshold, stop
        int min_in_row = *std::min_element(curr.begin(), curr.end());
        if (min_in_row > max_distance) {
            return max_distance + 1;
        }
        
        // Swap rows
        prev.swap(curr);
    }
    
    return prev[n];
}

