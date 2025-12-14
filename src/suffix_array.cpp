#include "suffix_array.h"
#include <cctype>

SuffixArray::SuffixArray() {
}

void SuffixArray::build(const std::string& text) {
    text_ = text;
    buildDoubling();
}

int SuffixArray::compareSuffix(int i, int j, int k) {
    if (rank_[i] != rank_[j]) {
        return rank_[i] - rank_[j];
    }
    
    int ri = (i + k < static_cast<int>(text_.length())) ? rank_[i + k] : -1;
    int rj = (j + k < static_cast<int>(text_.length())) ? rank_[j + k] : -1;
    
    return ri - rj;
}

void SuffixArray::buildDoubling() {
    int n = static_cast<int>(text_.length());
    
    // Initialize
    suffix_array_.resize(n);
    rank_.resize(n);
    
    for (int i = 0; i < n; ++i) {
        suffix_array_[i] = i;
        rank_[i] = std::toupper(text_[i]);
    }
    
    // Sort by first character
    std::sort(suffix_array_.begin(), suffix_array_.end(),
        [this](int a, int b) {
            return std::toupper(text_[a]) < std::toupper(text_[b]);
        });
    
    // Update ranks
    int r = 0;
    rank_[suffix_array_[0]] = 0;
    for (int i = 1; i < n; ++i) {
        if (std::toupper(text_[suffix_array_[i]]) != std::toupper(text_[suffix_array_[i-1]])) {
            r++;
        }
        rank_[suffix_array_[i]] = r;
    }
    
    // Doubling method
    std::vector<int> new_rank(n);
    for (int k = 1; k < n; k *= 2) {
        // Sort by (rank[i], rank[i+k])
        std::sort(suffix_array_.begin(), suffix_array_.end(),
            [this, k, n](int a, int b) {
                if (rank_[a] != rank_[b]) {
                    return rank_[a] < rank_[b];
                }
                int ra = (a + k < n) ? rank_[a + k] : -1;
                int rb = (b + k < n) ? rank_[b + k] : -1;
                return ra < rb;
            });
        
        // Update ranks
        new_rank[suffix_array_[0]] = 0;
        r = 0;
        for (int i = 1; i < n; ++i) {
            int a = suffix_array_[i-1];
            int b = suffix_array_[i];
            int ra = (a + k < n) ? rank_[a + k] : -1;
            int rb = (b + k < n) ? rank_[b + k] : -1;
            
            if (rank_[a] != rank_[b] || ra != rb) {
                r++;
            }
            new_rank[b] = r;
        }
        
        rank_ = new_rank;
        
        if (r == n - 1) break;  // All distinct
    }
}

std::pair<int, int> SuffixArray::binarySearch(const std::string& pattern) {
    int n = static_cast<int>(text_.length());
    int m = static_cast<int>(pattern.length());
    
    int left = 0, right = n - 1;
    int first = -1, last = -1;
    
    // Find first occurrence
    while (left <= right) {
        int mid = (left + right) / 2;
        int cmp = comparePatternSuffix(pattern, suffix_array_[mid]);
        
        if (cmp == 0) {
            first = mid;
            right = mid - 1;
        } else if (cmp < 0) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    
    if (first == -1) {
        return {-1, -1};
    }
    
    // Find last occurrence
    left = first;
    right = n - 1;
    while (left <= right) {
        int mid = (left + right) / 2;
        int cmp = comparePatternSuffix(pattern, suffix_array_[mid]);
        
        if (cmp == 0) {
            last = mid;
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    return {first, last};
}

int SuffixArray::comparePatternSuffix(const std::string& pattern, int suffix_idx) {
    int n = static_cast<int>(text_.length());
    int m = static_cast<int>(pattern.length());
    
    for (int i = 0; i < m; ++i) {
        if (suffix_idx + i >= n) {
            return 1;  // Pattern longer than suffix
        }
        
        char pc = std::toupper(pattern[i]);
        char sc = std::toupper(text_[suffix_idx + i]);
        
        if (pc < sc) return -1;
        if (pc > sc) return 1;
    }
    
    return 0;  // Match
}

SearchResult SuffixArray::search(const std::string& pattern) {
    SearchResult result;
    
    if (pattern.empty() || text_.empty()) {
        return result;
    }
    
    auto [first, last] = binarySearch(pattern);
    
    if (first == -1) {
        return result;
    }
    
    for (int i = first; i <= last; ++i) {
        result.positions.push_back(suffix_array_[i]);
    }
    
    result.count = static_cast<int>(result.positions.size());
    return result;
}

int SuffixArray::findFirst(const std::string& pattern) {
    auto [first, last] = binarySearch(pattern);
    if (first == -1) {
        return -1;
    }
    return suffix_array_[first];
}

std::vector<int> SuffixArray::buildLCP() {
    int n = static_cast<int>(text_.length());
    std::vector<int> lcp(n, 0);
    std::vector<int> inv_suffix(n);
    
    // Build inverse suffix array
    for (int i = 0; i < n; ++i) {
        inv_suffix[suffix_array_[i]] = i;
    }
    
    int k = 0;
    for (int i = 0; i < n; ++i) {
        if (inv_suffix[i] == n - 1) {
            k = 0;
            continue;
        }
        
        int j = suffix_array_[inv_suffix[i] + 1];
        
        while (i + k < n && j + k < n &&
               std::toupper(text_[i + k]) == std::toupper(text_[j + k])) {
            k++;
        }
        
        lcp[inv_suffix[i]] = k;
        
        if (k > 0) k--;
    }
    
    return lcp;
}

std::string SuffixArray::longestRepeatedSubstring() {
    std::vector<int> lcp = buildLCP();
    
    if (lcp.empty()) {
        return "";
    }
    
    int max_len = 0;
    int max_idx = 0;
    
    for (size_t i = 0; i < lcp.size(); ++i) {
        if (lcp[i] > max_len) {
            max_len = lcp[i];
            max_idx = i;
        }
    }
    
    if (max_len == 0) {
        return "";
    }
    
    return text_.substr(suffix_array_[max_idx], max_len);
}

