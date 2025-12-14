#include "bitap.h"
#include <algorithm>
#include <cctype>
#include <climits>

Bitap::Bitap(int max_errors) : max_errors_(max_errors) {
    if (max_errors_ < 0) max_errors_ = 2;
}

std::map<char, unsigned long long> Bitap::buildCharMask(const std::string& pattern) {
    std::map<char, unsigned long long> char_mask;
    
    unsigned long long mask = 1ULL;
    for (char c : pattern) {
        char upper_c = std::toupper(c);
        char_mask[upper_c] |= mask;
        mask <<= 1;
    }
    
    return char_mask;
}

std::vector<size_t> Bitap::shiftOrSearch(const std::string& text,
                                        const std::string& pattern,
                                        int max_errors) {
    std::vector<size_t> matches;
    
    if (text.empty() || pattern.empty() || pattern.length() > 64) {
        return matches;  // Bitap limited to 64 characters
    }
    
    int m = static_cast<int>(pattern.length());
    int n = static_cast<int>(text.length());
    
    if (max_errors < 0) {
        max_errors = max_errors_;
    }
    
    // Build character mask
    std::map<char, unsigned long long> char_mask = buildCharMask(pattern);
    
    // Initialize state vectors for each error level
    std::vector<unsigned long long> R(max_errors + 1, ~0ULL);
    
    // Pattern mask (all bits set except last)
    unsigned long long pattern_mask = ~(1ULL << (m - 1));
    
    for (int i = 0; i < n; ++i) {
        char c = std::toupper(text[i]);
        unsigned long long char_mask_val = 0;
        
        if (char_mask.find(c) != char_mask.end()) {
            char_mask_val = char_mask[c];
        }
        
        // Update state for each error level
        unsigned long long old_R0 = R[0];
        
        // Exact match (0 errors)
        R[0] = (R[0] << 1) | char_mask_val;
        R[0] |= pattern_mask;
        
        // With errors
        for (int e = 1; e <= max_errors; ++e) {
            unsigned long long new_R = (R[e] << 1) | char_mask_val;
            
            // Substitution
            new_R |= (old_R0 << 1);
            
            // Insertion
            new_R |= R[e-1];
            
            // Deletion
            new_R |= (R[e-1] << 1);
            
            old_R0 = R[e];
            R[e] = new_R | pattern_mask;
        }
        
        // Check for match (bit 0 of R[max_errors] is 0)
        if ((R[max_errors] & 1) == 0) {
            matches.push_back(i - m + 1);
        }
    }
    
    return matches;
}

SearchResult Bitap::search(const std::string& text,
                           const std::string& pattern,
                           int max_errors) {
    SearchResult result;
    
    if (max_errors < 0) {
        max_errors = max_errors_;
    }
    
    std::vector<size_t> positions = shiftOrSearch(text, pattern, max_errors);
    result.positions.clear();
    for (size_t pos : positions) {
        result.positions.push_back(static_cast<int>(pos));
    }
    result.count = static_cast<int>(positions.size());
    
    return result;
}

bool Bitap::contains(const std::string& text,
                    const std::string& pattern,
                    int max_errors) {
    return findFirst(text, pattern, max_errors) != -1;
}

int Bitap::findFirst(const std::string& text,
                    const std::string& pattern,
                    int max_errors) {
    std::vector<size_t> positions = shiftOrSearch(text, pattern, max_errors);
    if (positions.empty()) {
        return -1;
    }
    return static_cast<int>(positions[0]);
}

int Bitap::bitapEditDistance(const std::string& text, const std::string& pattern) {
    // Simplified edit distance using bit-parallelism
    // This is a basic implementation
    
    if (text.empty()) return static_cast<int>(pattern.length());
    if (pattern.empty()) return static_cast<int>(text.length());
    
    // For longer sequences, use standard edit distance
    if (pattern.length() > 64) {
        // Fallback to standard method
        int m = static_cast<int>(text.length());
        int n = static_cast<int>(pattern.length());
        std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));
        
        for (int i = 0; i <= m; ++i) dp[i][0] = i;
        for (int j = 0; j <= n; ++j) dp[0][j] = j;
        
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                int cost = (std::toupper(text[i-1]) == std::toupper(pattern[j-1])) ? 0 : 1;
                dp[i][j] = std::min({dp[i-1][j] + 1, dp[i][j-1] + 1, dp[i-1][j-1] + cost});
            }
        }
        
        return dp[m][n];
    }
    
    // Bit-parallel approach for short patterns
    std::map<char, unsigned long long> char_mask = buildCharMask(pattern);
    int m = static_cast<int>(pattern.length());
    int n = static_cast<int>(text.length());
    
    unsigned long long R = ~0ULL;
    unsigned long long pattern_mask = ~(1ULL << (m - 1));
    
    for (int i = 0; i < n; ++i) {
        char c = std::toupper(text[i]);
        unsigned long long char_mask_val = 0;
        
        if (char_mask.find(c) != char_mask.end()) {
            char_mask_val = char_mask[c];
        }
        
        R = ((R << 1) | char_mask_val) | pattern_mask;
    }
    
    // Count errors (simplified)
    int errors = 0;
    unsigned long long temp = R;
    while (temp != 0) {
        if (temp & 1) errors++;
        temp >>= 1;
    }
    
    return errors;
}

