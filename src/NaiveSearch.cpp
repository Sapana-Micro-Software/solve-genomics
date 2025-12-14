#include "NaiveSearch.h"
#include "utils.h"
#include <algorithm>

SearchResult NaiveSearch::search(const std::string& sequence, const std::string& pattern) {
    SearchResult result;                                                                      // Initialize result structure
    
    if (pattern.empty() || sequence.empty() || pattern.length() > sequence.length()) {        // Validate input parameters
        return result;                                                                        // Return empty result if invalid
    }
    
    size_t n = sequence.length();                                                             // Get sequence length
    size_t m = pattern.length();                                                              // Get pattern length
    
    // Slide pattern over sequence one position at a time
    for (size_t i = 0; i <= n - m; ++i) {                                                     // Slide window through sequence
        size_t j = 0;                                                                         // Start comparing from beginning of pattern
        
        // Compare pattern with current window
        while (j < m && std::toupper(sequence[i + j]) == std::toupper(pattern[j])) {          // Continue while characters match
            j++;                                                                              // Advance to next character
        }
        
        // If we matched all characters, pattern found
        if (j == m) {                                                                         // Check if entire pattern matched
            result.positions.push_back(static_cast<int>(i));                                 // Record position of match
            result.count++;                                                                   // Increment match counter
        }
    }
    
    return result;                                                                            // Return all found positions
}

int NaiveSearch::findFirst(const std::string& sequence, const std::string& pattern) {
    if (pattern.empty() || sequence.empty() || pattern.length() > sequence.length()) {        // Validate input parameters
        return -1;                                                                            // Return -1 if invalid
    }
    
    size_t n = sequence.length();                                                             // Get sequence length
    size_t m = pattern.length();                                                              // Get pattern length
    
    for (size_t i = 0; i <= n - m; ++i) {                                                     // Slide window through sequence
        size_t j = 0;                                                                         // Start comparing from beginning of pattern
        while (j < m && std::toupper(sequence[i + j]) == std::toupper(pattern[j])) {          // Continue while characters match
            j++;                                                                              // Advance to next character
        }
        
        if (j == m) {                                                                         // Check if entire pattern matched
            return static_cast<int>(i);                                                       // Return first matching position
        }
    }
    
    return -1;                                                                                // Pattern not found
}

