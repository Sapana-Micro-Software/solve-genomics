#include "ExactMatch.h"
#include "utils.h"
#include <algorithm>

SearchResult ExactMatch::search(const std::string& sequence, const std::string& pattern) {
    SearchResult result;                                                                      // Initialize result structure
    
    if (pattern.empty() || sequence.empty() || pattern.length() > sequence.length()) {        // Validate input parameters
        return result;                                                                        // Return empty result if invalid
    }
    
    size_t seq_len = sequence.length();                                                       // Get sequence length
    size_t pat_len = pattern.length();                                                        // Get pattern length
    
    // Linear search through the sequence
    for (size_t i = 0; i <= seq_len - pat_len; ++i) {                                         // Slide window through sequence
        bool match = true;                                                                    // Assume match until proven otherwise
        for (size_t j = 0; j < pat_len; ++j) {                                                // Compare each character in pattern
            // Case-insensitive comparison
            if (std::toupper(sequence[i + j]) != std::toupper(pattern[j])) {                  // Compare characters (case-insensitive)
                match = false;                                                                // Mismatch found, stop checking
                break;
            }
        }
        
        if (match) {                                                                          // If all characters matched
            result.positions.push_back(static_cast<int>(i));                                 // Record position of match
            result.count++;                                                                   // Increment match counter
        }
    }
    
    return result;                                                                            // Return all found positions
}

bool ExactMatch::contains(const std::string& sequence, const std::string& pattern) {
    return findFirst(sequence, pattern) != -1;                                                // Check if pattern exists in sequence
}

int ExactMatch::findFirst(const std::string& sequence, const std::string& pattern) {
    if (pattern.empty() || sequence.empty() || pattern.length() > sequence.length()) {        // Validate input parameters
        return -1;                                                                            // Return -1 if invalid
    }
    
    size_t seq_len = sequence.length();                                                       // Get sequence length
    size_t pat_len = pattern.length();                                                        // Get pattern length
    
    for (size_t i = 0; i <= seq_len - pat_len; ++i) {                                         // Slide window through sequence
        bool match = true;                                                                    // Assume match until proven otherwise
        for (size_t j = 0; j < pat_len; ++j) {                                                // Compare each character in pattern
            if (std::toupper(sequence[i + j]) != std::toupper(pattern[j])) {                  // Compare characters (case-insensitive)
                match = false;                                                                // Mismatch found, stop checking
                break;
            }
        }
        
        if (match) {                                                                          // If all characters matched
            return static_cast<int>(i);                                                       // Return first matching position
        }
    }
    
    return -1;                                                                                // Pattern not found
}

