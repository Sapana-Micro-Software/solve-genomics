#include "classic_string_matching.h"
#include <algorithm>
#include <cmath>
#include <cctype>

// Rabin-Karp Implementation

SearchResult ClassicStringMatching::RabinKarp::search(const std::string& text, const std::string& pattern) {
    return searchWithHash(text, pattern, 256, 101);
}

SearchResult ClassicStringMatching::RabinKarp::searchWithHash(const std::string& text,
                                                             const std::string& pattern,
                                                             int base,
                                                             int modulus) {
    SearchResult result;                                                                      // Initialize result structure
    
    if (pattern.empty() || text.empty() || pattern.length() > text.length()) {                // Validate input parameters
        return result;                                                                        // Return empty result if invalid
    }
    
    int pattern_len = pattern.length();                                                       // Get pattern length
    int text_len = text.length();                                                             // Get text length
    
    // Calculate hash for pattern and first window of text
    int pattern_hash = calculateHash(pattern, base, modulus);                                 // Compute hash value for pattern
    int text_hash = calculateHash(text.substr(0, pattern_len), base, modulus);               // Compute hash for first window
    
    // Calculate base^(pattern_len-1) for rolling hash
    int base_power = 1;                                                                       // Initialize power calculation
    for (int i = 0; i < pattern_len - 1; ++i) {                                               // Calculate base^(pattern_len-1)
        base_power = (base_power * base) % modulus;                                           // Modular exponentiation
    }
    
    // Slide pattern over text
    for (int i = 0; i <= text_len - pattern_len; ++i) {                                       // Slide window through text
        // Check hash values
        if (pattern_hash == text_hash) {                                                      // Compare hash values (fast check)
            // Hash match - verify with character-by-character comparison
            bool match = true;                                                                // Assume match until proven otherwise
            for (int j = 0; j < pattern_len; ++j) {                                           // Verify character-by-character
                if (std::toupper(text[i + j]) != std::toupper(pattern[j])) {                  // Compare characters (case-insensitive)
                    match = false;                                                            // Mismatch found, stop checking
                    break;
                }
            }
            
            if (match) {                                                                      // If all characters matched
                result.positions.push_back(i);                                               // Record position of match
                result.count++;                                                               // Increment match counter
            }
        }
        
        // Calculate hash for next window (rolling hash)
        if (i < text_len - pattern_len) {                                                     // Check if more windows remain
            text_hash = rollingHash(text_hash, text[i], text[i + pattern_len],
                                   base, modulus, pattern_len);                               // Update hash for next window
            // Adjust for negative modulo
            if (text_hash < 0) {                                                              // Handle negative modulo result
                text_hash += modulus;                                                         // Make positive
            }
        }
    }
    
    return result;                                                                            // Return all found positions
}

int ClassicStringMatching::RabinKarp::calculateHash(const std::string& str, int base, int modulus) {
    int hash = 0;                                                                             // Initialize hash value
    for (char c : str) {                                                                      // Iterate through each character
        hash = (hash * base + std::toupper(c)) % modulus;                                    // Polynomial rolling hash calculation
    }
    return hash;                                                                              // Return computed hash value
}

int ClassicStringMatching::RabinKarp::rollingHash(int old_hash, char old_char, char new_char,
                                                  int base, int modulus, int pattern_len) {
    // Remove old character and add new character
    int base_power = 1;                                                                       // Initialize power calculation
    for (int i = 0; i < pattern_len - 1; ++i) {                                               // Calculate base^(pattern_len-1)
        base_power = (base_power * base) % modulus;                                           // Modular exponentiation
    }
    
    int new_hash = (old_hash - std::toupper(old_char) * base_power) % modulus;               // Remove contribution of old character
    new_hash = (new_hash * base + std::toupper(new_char)) % modulus;                         // Add contribution of new character
    
    return new_hash;                                                                          // Return updated hash value
}

// KMP Implementation

SearchResult ClassicStringMatching::KMP::search(const std::string& text, const std::string& pattern) {
    SearchResult result;                                                                      // Initialize result structure
    
    if (pattern.empty() || text.empty() || pattern.length() > text.length()) {                // Validate input parameters
        return result;                                                                        // Return empty result if invalid
    }
    
    std::vector<int> lps = buildFailureFunction(pattern);                                     // Build longest proper prefix-suffix array
    
    int text_len = text.length();                                                             // Get text length
    int pattern_len = pattern.length();                                                       // Get pattern length
    int i = 0;                                                                                // Index for text (never decreases)
    int j = 0;                                                                                // Index for pattern (can decrease via lps)
    
    while (i < text_len) {                                                                    // Continue until end of text
        if (std::toupper(text[i]) == std::toupper(pattern[j])) {                              // Characters match
            i++;                                                                              // Advance text index
            j++;                                                                              // Advance pattern index
        }
        
        if (j == pattern_len) {                                                               // Entire pattern matched
            // Pattern found
            result.positions.push_back(i - j);                                               // Record position of match
            result.count++;                                                                   // Increment match counter
            j = lps[j - 1];                                                                   // Continue searching using failure function
        } else if (i < text_len && std::toupper(text[i]) != std::toupper(pattern[j])) {       // Mismatch found
            if (j != 0) {                                                                     // If pattern index not at start
                j = lps[j - 1];                                                               // Use failure function to skip characters
            } else {
                i++;                                                                          // Advance text index if at pattern start
            }
        }
    }
    
    return result;                                                                            // Return all found positions
}

int ClassicStringMatching::KMP::findFirst(const std::string& text, const std::string& pattern) {
    if (pattern.empty() || text.empty() || pattern.length() > text.length()) {
        return -1;
    }
    
    std::vector<int> lps = buildFailureFunction(pattern);
    
    int text_len = text.length();
    int pattern_len = pattern.length();
    int i = 0;
    int j = 0;
    
    while (i < text_len) {
        if (std::toupper(text[i]) == std::toupper(pattern[j])) {
            i++;
            j++;
        }
        
        if (j == pattern_len) {
            return i - j;  // Found first match
        } else if (i < text_len && std::toupper(text[i]) != std::toupper(pattern[j])) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
    
    return -1;
}

std::vector<int> ClassicStringMatching::KMP::buildFailureFunction(const std::string& pattern) {
    int pattern_len = pattern.length();                                                        // Get pattern length
    std::vector<int> lps(pattern_len, 0);                                                     // Initialize LPS array (Longest Proper Prefix-Suffix)
    
    int len = 0;                                                                              // Length of previous longest prefix suffix
    int i = 1;                                                                                // Start from second character
    
    while (i < pattern_len) {                                                                 // Build LPS array for each position
        if (std::toupper(pattern[i]) == std::toupper(pattern[len])) {                         // Characters match
            len++;                                                                            // Extend matching prefix-suffix
            lps[i] = len;                                                                     // Store length at current position
            i++;                                                                              // Move to next character
        } else {
            if (len != 0) {                                                                   // If previous match exists
                len = lps[len - 1];                                                           // Use previous LPS value (don't increment i)
            } else {
                lps[i] = 0;                                                                   // No matching prefix-suffix
                i++;                                                                          // Move to next character
            }
        }
    }
    
    return lps;                                                                               // Return failure function array
}

// Boyer-Moore Implementation

SearchResult ClassicStringMatching::BoyerMoore::search(const std::string& text, const std::string& pattern) {
    SearchResult result;                                                                      // Initialize result structure
    
    if (pattern.empty() || text.empty() || pattern.length() > text.length()) {                // Validate input parameters
        return result;                                                                        // Return empty result if invalid
    }
    
    std::vector<int> bad_char = buildBadCharacterTable(pattern);                             // Build bad character shift table
    std::vector<int> good_suffix = buildGoodSuffixTable(pattern);                            // Build good suffix shift table
    
    int text_len = text.length();                                                             // Get text length
    int pattern_len = pattern.length();                                                       // Get pattern length
    int shift = 0;                                                                            // Current shift position in text
    
    while (shift <= text_len - pattern_len) {                                                 // Continue while pattern can fit
        int j = pattern_len - 1;                                                              // Start matching from end of pattern
        
        // Match from right to left
        while (j >= 0 && std::toupper(text[shift + j]) == std::toupper(pattern[j])) {        // Compare characters right-to-left
            j--;                                                                              // Continue matching backwards
        }
        
        if (j < 0) {                                                                          // Entire pattern matched
            // Pattern found
            result.positions.push_back(shift);                                               // Record position of match
            result.count++;                                                                   // Increment match counter
            
            // Shift by good suffix
            shift += (shift + pattern_len < text_len) ?                                      // Check if more text remains
                    good_suffix[0] : 1;                                                       // Use good suffix shift or shift by 1
        } else {
            // Bad character heuristic
            char bad_char_val = std::toupper(text[shift + j]);                               // Get mismatched character from text
            int bad_char_shift = 1;                                                           // Default shift of 1
            if (bad_char_val >= 0 && bad_char_val < 128) {                                   // Check if valid ASCII character
                int char_pos = bad_char[bad_char_val];                                        // Get rightmost position of character in pattern
                if (char_pos < pattern_len) {                                                 // If character exists in pattern
                    bad_char_shift = std::max(1, j - char_pos);                              // Calculate shift to align character
                }
            }
            
            // Good suffix heuristic
            int good_suffix_shift = good_suffix[j + 1];                                      // Get shift from good suffix table
            
            // Take maximum shift
            shift += std::max(bad_char_shift, good_suffix_shift);                            // Use larger shift for efficiency
        }
    }
    
    return result;                                                                            // Return all found positions
}

int ClassicStringMatching::BoyerMoore::findFirst(const std::string& text, const std::string& pattern) {
    if (pattern.empty() || text.empty() || pattern.length() > text.length()) {
        return -1;
    }
    
    std::vector<int> bad_char = buildBadCharacterTable(pattern);
    std::vector<int> good_suffix = buildGoodSuffixTable(pattern);
    
    int text_len = text.length();
    int pattern_len = pattern.length();
    int shift = 0;
    
    while (shift <= text_len - pattern_len) {
        int j = pattern_len - 1;
        
        while (j >= 0 && std::toupper(text[shift + j]) == std::toupper(pattern[j])) {
            j--;
        }
        
        if (j < 0) {
            return shift;  // Found first match
        } else {
            char bad_char_val = std::toupper(text[shift + j]);
            int bad_char_shift = 1;
            if (bad_char_val >= 0 && bad_char_val < 128) {
                int char_pos = bad_char[bad_char_val];
                if (char_pos < pattern_len) {
                    bad_char_shift = std::max(1, j - char_pos);
                }
            }
            int good_suffix_shift = good_suffix[j + 1];
            
            shift += std::max(bad_char_shift, good_suffix_shift);
        }
    }
    
    return -1;
}

std::vector<int> ClassicStringMatching::BoyerMoore::buildBadCharacterTable(const std::string& pattern) {
    const int ALPHABET_SIZE = 128;  // ASCII
    std::vector<int> bad_char(ALPHABET_SIZE, -1);
    
    // Fill with rightmost occurrence of each character
    for (int i = 0; i < static_cast<int>(pattern.length()); ++i) {
        char c = std::toupper(pattern[i]);
        if (c >= 0 && c < ALPHABET_SIZE) {
            bad_char[c] = i;
        }
    }
    
    return bad_char;
}

std::vector<int> ClassicStringMatching::BoyerMoore::buildGoodSuffixTable(const std::string& pattern) {
    int pattern_len = pattern.length();
    std::vector<int> good_suffix(pattern_len + 1, 0);
    std::vector<int> suffix = buildSuffixArray(pattern);
    
    // Case 1: Pattern matches after shift
    for (int i = 0; i < pattern_len; ++i) {
        good_suffix[i] = pattern_len;
    }
    
    // Case 2: Suffix of pattern matches prefix
    int j = 0;
    for (int i = pattern_len - 1; i >= 0; --i) {
        if (suffix[i] == i + 1) {
            for (; j < pattern_len - 1 - i; ++j) {
                if (good_suffix[j] == pattern_len) {
                    good_suffix[j] = pattern_len - 1 - i;
                }
            }
        }
    }
    
    // Case 3: Suffix of pattern matches suffix of text
    for (int i = 0; i < pattern_len - 1; ++i) {
        int suffix_len = suffix[i];
        good_suffix[pattern_len - 1 - suffix_len] = pattern_len - 1 - i;
    }
    
    return good_suffix;
}

std::vector<int> ClassicStringMatching::BoyerMoore::buildSuffixArray(const std::string& pattern) {
    int pattern_len = pattern.length();
    std::vector<int> suffix(pattern_len, 0);
    
    suffix[pattern_len - 1] = pattern_len;
    int g = pattern_len - 1;
    int f = 0;
    
    for (int i = pattern_len - 2; i >= 0; --i) {
        if (i > g && suffix[i + pattern_len - 1 - f] < i - g) {
            suffix[i] = suffix[i + pattern_len - 1 - f];
        } else {
            if (i < g) {
                g = i;
            }
            f = i;
            
            while (g >= 0 && std::toupper(pattern[g]) == 
                   std::toupper(pattern[g + pattern_len - 1 - f])) {
                g--;
            }
            suffix[i] = f - g;
        }
    }
    
    return suffix;
}

// Comparison function

ClassicStringMatching::ComparisonResult ClassicStringMatching::compareAll(const std::string& text, 
                                                                          const std::string& pattern) {
    ComparisonResult result;
    result.rabin_karp = RabinKarp::search(text, pattern);
    result.kmp = KMP::search(text, pattern);
    result.boyer_moore = BoyerMoore::search(text, pattern);
    return result;
}

