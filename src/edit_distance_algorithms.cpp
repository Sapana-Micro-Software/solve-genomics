#include "edit_distance_algorithms.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <cctype>

int EditDistanceAlgorithms::levenshteinDistance(const std::string& str1, const std::string& str2) {
    int m = str1.length();
    int n = str2.length();
    
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    // Initialize
    for (int i = 0; i <= m; ++i) dp[i][0] = i;
    for (int j = 0; j <= n; ++j) dp[0][j] = j;
    
    // Fill matrix
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + std::min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
            }
        }
    }
    
    return dp[m][n];
}

int EditDistanceAlgorithms::levenshteinDistanceOptimized(const std::string& str1, const std::string& str2) {
    int m = str1.length();
    int n = str2.length();
    
    // Use shorter string for columns to minimize space
    if (m < n) {
        return levenshteinDistanceOptimized(str2, str1);
    }
    
    std::vector<int> prev(n + 1);
    std::vector<int> curr(n + 1);
    
    for (int j = 0; j <= n; ++j) {
        prev[j] = j;
    }
    
    for (int i = 1; i <= m; ++i) {
        curr[0] = i;
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                curr[j] = prev[j - 1];
            } else {
                curr[j] = 1 + std::min({prev[j], curr[j - 1], prev[j - 1]});
            }
        }
        prev.swap(curr);
    }
    
    return prev[n];
}

int EditDistanceAlgorithms::levenshteinDistanceBounded(const std::string& str1,
                                                      const std::string& str2,
                                                      int max_distance) {
    int m = str1.length();
    int n = str2.length();
    
    // Early termination
    int len_diff = (m > n) ? (m - n) : (n - m);
    if (len_diff > max_distance) {
        return max_distance + 1;
    }
    
    // Use optimized version with early termination
    std::vector<int> prev(n + 1);
    std::vector<int> curr(n + 1);
    
    for (int j = 0; j <= n; ++j) {
        prev[j] = j;
    }
    
    for (int i = 1; i <= m; ++i) {
        curr[0] = i;
        bool found_valid = false;
        
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                curr[j] = prev[j - 1];
            } else {
                curr[j] = 1 + std::min({prev[j], curr[j - 1], prev[j - 1]});
            }
            
            if (curr[j] <= max_distance) {
                found_valid = true;
            }
        }
        
        if (!found_valid) {
            return max_distance + 1;
        }
        
        prev.swap(curr);
    }
    
    return prev[n];
}

int EditDistanceAlgorithms::damerauLevenshteinDistance(const std::string& str1, const std::string& str2) {
    int m = str1.length();
    int n = str2.length();
    
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    for (int i = 0; i <= m; ++i) dp[i][0] = i;
    for (int j = 0; j <= n; ++j) dp[0][j] = j;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + std::min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
            }
            
            // Check for transposition
            if (i > 1 && j > 1 &&
                std::toupper(str1[i - 1]) == std::toupper(str2[j - 2]) &&
                std::toupper(str1[i - 2]) == std::toupper(str2[j - 1])) {
                dp[i][j] = std::min(dp[i][j], dp[i - 2][j - 2] + 1);
            }
        }
    }
    
    return dp[m][n];
}

int EditDistanceAlgorithms::hammingDistance(const std::string& str1, const std::string& str2) {
    if (str1.length() != str2.length()) {
        return -1;
    }
    
    int distance = 0;
    for (size_t i = 0; i < str1.length(); ++i) {
        if (std::toupper(str1[i]) != std::toupper(str2[i])) {
            distance++;
        }
    }
    
    return distance;
}

double EditDistanceAlgorithms::jaroWinklerDistance(const std::string& str1, const std::string& str2) {
    if (str1 == str2) return 1.0;
    if (str1.empty() || str2.empty()) return 0.0;
    
    // Jaro distance
    int match_window = std::max(str1.length(), str2.length()) / 2 - 1;
    if (match_window < 0) match_window = 0;
    
    std::vector<bool> str1_matches(str1.length(), false);
    std::vector<bool> str2_matches(str2.length(), false);
    
    int matches = 0;
    int transpositions = 0;
    
    // Find matches
    for (size_t i = 0; i < str1.length(); ++i) {
        int start = (i >= static_cast<size_t>(match_window)) ? 
                    (i - match_window) : 0;
        int end = std::min(i + match_window + 1, str2.length());
        
        for (int j = start; j < end; ++j) {
            if (str2_matches[j] || std::toupper(str1[i]) != std::toupper(str2[j])) {
                continue;
            }
            str1_matches[i] = true;
            str2_matches[j] = true;
            matches++;
            break;
        }
    }
    
    if (matches == 0) return 0.0;
    
    // Count transpositions
    int k = 0;
    for (size_t i = 0; i < str1.length(); ++i) {
        if (!str1_matches[i]) continue;
        while (!str2_matches[k]) k++;
        if (std::toupper(str1[i]) != std::toupper(str2[k])) {
            transpositions++;
        }
        k++;
    }
    
    double jaro = (matches / static_cast<double>(str1.length()) +
                   matches / static_cast<double>(str2.length()) +
                   (matches - transpositions / 2.0) / matches) / 3.0;
    
    // Winkler modification (common prefix)
    int prefix_len = 0;
    int max_prefix = std::min(4, static_cast<int>(std::min(str1.length(), str2.length())));
    for (int i = 0; i < max_prefix; ++i) {
        if (std::toupper(str1[i]) == std::toupper(str2[i])) {
            prefix_len++;
        } else {
            break;
        }
    }
    
    return jaro + (0.1 * prefix_len * (1 - jaro));
}

int EditDistanceAlgorithms::longestCommonSubsequence(const std::string& str1, const std::string& str2) {
    int m = str1.length();
    int n = str2.length();
    
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
                dp[i][j] = std::max(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }
    
    return dp[m][n];
}

int EditDistanceAlgorithms::longestCommonSubstring(const std::string& str1, const std::string& str2) {
    int m = str1.length();
    int n = str2.length();
    int max_len = 0;
    
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                max_len = std::max(max_len, dp[i][j]);
            } else {
                dp[i][j] = 0;
            }
        }
    }
    
    return max_len;
}

int EditDistanceAlgorithms::weightedEditDistance(const std::string& str1,
                                                 const std::string& str2,
                                                 int insert_cost,
                                                 int delete_cost,
                                                 int substitute_cost) {
    int m = str1.length();
    int n = str2.length();
    
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    for (int i = 0; i <= m; ++i) dp[i][0] = i * delete_cost;
    for (int j = 0; j <= n; ++j) dp[0][j] = j * insert_cost;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = std::min({
                    dp[i - 1][j] + delete_cost,
                    dp[i][j - 1] + insert_cost,
                    dp[i - 1][j - 1] + substitute_cost
                });
            }
        }
    }
    
    return dp[m][n];
}

bool EditDistanceAlgorithms::isTransition(char base1, char base2) {
    char b1 = std::toupper(base1);
    char b2 = std::toupper(base2);
    
    // Transitions: A<->G (purines) or T<->C (pyrimidines)
    return (b1 == 'A' && b2 == 'G') || (b1 == 'G' && b2 == 'A') ||
           (b1 == 'T' && b2 == 'C') || (b1 == 'C' && b2 == 'T');
}

int EditDistanceAlgorithms::dnaEditDistance(const std::string& seq1, const std::string& seq2) {
    // DNA-specific: transitions (A<->G, T<->C) cost less than transversions
    int transition_cost = 1;
    int transversion_cost = 2;
    int insert_cost = 2;
    int delete_cost = 2;
    
    int m = seq1.length();
    int n = seq2.length();
    
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    for (int i = 0; i <= m; ++i) dp[i][0] = i * delete_cost;
    for (int j = 0; j <= n; ++j) dp[0][j] = j * insert_cost;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            char b1 = std::toupper(seq1[i - 1]);
            char b2 = std::toupper(seq2[j - 1]);
            
            if (b1 == b2) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                int sub_cost = isTransition(b1, b2) ? transition_cost : transversion_cost;
                dp[i][j] = std::min({
                    dp[i - 1][j] + delete_cost,
                    dp[i][j - 1] + insert_cost,
                    dp[i - 1][j - 1] + sub_cost
                });
            }
        }
    }
    
    return dp[m][n];
}

