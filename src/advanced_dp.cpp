#include "advanced_dp.h"
#include <algorithm>
#include <cctype>
#include <climits>

AlignmentResult AdvancedDP::smithWatermanSpaceOptimized(const std::string& seq1,
                                                        const std::string& seq2,
                                                        int match_score,
                                                        int mismatch_score,
                                                        int gap_penalty) {
    AlignmentResult result;
    
    if (seq1.empty() || seq2.empty()) {
        return result;
    }
    
    int m = static_cast<int>(seq1.length());
    int n = static_cast<int>(seq2.length());
    
    // Use shorter sequence for columns
    bool swap = (m > n);
    const std::string& s1 = swap ? seq2 : seq1;
    const std::string& s2 = swap ? seq1 : seq2;
    if (swap) {
        int temp = m;
        m = n;
        n = temp;
    }
    
    // Only two rows needed
    std::vector<int> prev_row(n + 1, 0);
    std::vector<int> curr_row(n + 1, 0);
    
    int max_score = 0;
    int max_i = 0, max_j = 0;
    
    for (int i = 1; i <= m; ++i) {
        curr_row[0] = 0;
        for (int j = 1; j <= n; ++j) {
            int match = prev_row[j-1] + score(s1[i-1], s2[j-1], match_score, mismatch_score);
            int del = prev_row[j] + gap_penalty;
            int ins = curr_row[j-1] + gap_penalty;
            
            curr_row[j] = std::max({0, match, del, ins});
            
            if (curr_row[j] > max_score) {
                max_score = curr_row[j];
                max_i = i;
                max_j = j;
            }
        }
        prev_row = curr_row;
    }
    
    result.score = max_score;
    
    // Simplified traceback (full traceback would require storing more information)
    // For space-optimized version, we'd need to recompute or use divide-and-conquer
    
    return result;
}

AlignmentResult AdvancedDP::needlemanWunschSpaceOptimized(const std::string& seq1,
                                                          const std::string& seq2,
                                                          int match_score,
                                                          int mismatch_score,
                                                          int gap_penalty) {
    AlignmentResult result;
    
    if (seq1.empty() || seq2.empty()) {
        return result;
    }
    
    int m = static_cast<int>(seq1.length());
    int n = static_cast<int>(seq2.length());
    
    // Two rows
    std::vector<int> prev_row(n + 1);
    std::vector<int> curr_row(n + 1);
    
    // Initialize first row
    for (int j = 0; j <= n; ++j) {
        prev_row[j] = j * gap_penalty;
    }
    
    for (int i = 1; i <= m; ++i) {
        curr_row[0] = i * gap_penalty;
        for (int j = 1; j <= n; ++j) {
            int match = prev_row[j-1] + score(seq1[i-1], seq2[j-1], match_score, mismatch_score);
            int del = prev_row[j] + gap_penalty;
            int ins = curr_row[j-1] + gap_penalty;
            
            curr_row[j] = std::max({match, del, ins});
        }
        prev_row = curr_row;
    }
    
    result.score = curr_row[n];
    
    return result;
}

AlignmentResult AdvancedDP::affineGapAlignment(const std::string& seq1,
                                              const std::string& seq2,
                                              int match_score,
                                              int mismatch_score,
                                              int gap_open,
                                              int gap_extend) {
    AlignmentResult result;
    
    if (seq1.empty() || seq2.empty()) {
        return result;
    }
    
    int m = static_cast<int>(seq1.length());
    int n = static_cast<int>(seq2.length());
    
    // Three matrices: M (match), I (insert), D (delete)
    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, INT_MIN));
    std::vector<std::vector<int>> I(m + 1, std::vector<int>(n + 1, INT_MIN));
    std::vector<std::vector<int>> D(m + 1, std::vector<int>(n + 1, INT_MIN));
    
    // Initialize
    M[0][0] = 0;
    for (int i = 1; i <= m; ++i) {
        I[i][0] = gap_open + (i - 1) * gap_extend;
    }
    for (int j = 1; j <= n; ++j) {
        D[0][j] = gap_open + (j - 1) * gap_extend;
    }
    
    // Fill matrices
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int s = score(seq1[i-1], seq2[j-1], match_score, mismatch_score);
            
            M[i][j] = s + std::max({M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]});
            I[i][j] = std::max({M[i][j-1] + gap_open, I[i][j-1] + gap_extend});
            D[i][j] = std::max({M[i-1][j] + gap_open, D[i-1][j] + gap_extend});
        }
    }
    
    result.score = std::max({M[m][n], I[m][n], D[m][n]});
    
    return result;
}

AlignmentResult AdvancedDP::bandedAlignment(const std::string& seq1,
                                            const std::string& seq2,
                                            int bandwidth,
                                            int match_score,
                                            int mismatch_score,
                                            int gap_penalty) {
    AlignmentResult result;
    
    if (seq1.empty() || seq2.empty()) {
        return result;
    }
    
    int m = static_cast<int>(seq1.length());
    int n = static_cast<int>(seq2.length());
    
    // Only compute within band
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, INT_MIN));
    
    dp[0][0] = 0;
    
    for (int i = 0; i <= m; ++i) {
        int start_j = std::max(0, i - bandwidth);
        int end_j = std::min(n, i + bandwidth);
        
        for (int j = start_j; j <= end_j; ++j) {
            if (i == 0 && j == 0) continue;
            
            int best = INT_MIN;
            
            if (i > 0 && j > 0) {
                int s = score(seq1[i-1], seq2[j-1], match_score, mismatch_score);
                best = std::max(best, dp[i-1][j-1] + s);
            }
            
            if (i > 0 && j >= start_j) {
                best = std::max(best, dp[i-1][j] + gap_penalty);
            }
            
            if (j > 0 && j <= end_j) {
                best = std::max(best, dp[i][j-1] + gap_penalty);
            }
            
            dp[i][j] = best;
        }
    }
    
    result.score = dp[m][n];
    
    return result;
}

AlignmentResult AdvancedDP::hirschbergAlignment(const std::string& seq1,
                                                const std::string& seq2,
                                                int match_score,
                                                int mismatch_score,
                                                int gap_penalty) {
    AlignmentResult result;
    
    if (seq1.empty() || seq2.empty()) {
        if (seq1.empty() && seq2.empty()) {
            result.score = 0;
        } else {
            result.score = static_cast<int>(std::max(seq1.length(), seq2.length())) * gap_penalty;
        }
        return result;
    }
    
    if (seq1.length() <= 1 || seq2.length() <= 1) {
        // Base case: use standard DP
        return needlemanWunschSpaceOptimized(seq1, seq2, match_score, mismatch_score, gap_penalty);
    }
    
    int mid = static_cast<int>(seq1.length()) / 2;
    
    // Find optimal midpoint
    int best_j = findMidpoint(seq1, seq2, match_score, mismatch_score, gap_penalty);
    
    // Recurse on left and right halves
    std::string left1 = seq1.substr(0, mid);
    std::string right1 = seq1.substr(mid);
    std::string left2 = seq2.substr(0, best_j);
    std::string right2 = seq2.substr(best_j);
    
    AlignmentResult left_result = hirschbergAlignment(left1, left2, match_score, mismatch_score, gap_penalty);
    AlignmentResult right_result = hirschbergAlignment(right1, right2, match_score, mismatch_score, gap_penalty);
    
    result.score = left_result.score + right_result.score;
    
    return result;
}

int AdvancedDP::score(char a, char b, int match_score, int mismatch_score) {
    if (std::toupper(a) == std::toupper(b)) {
        return match_score;
    }
    return mismatch_score;
}

int AdvancedDP::findMidpoint(const std::string& seq1, const std::string& seq2,
                             int match_score, int mismatch_score, int gap_penalty) {
    int mid = static_cast<int>(seq1.length()) / 2;
    int n = static_cast<int>(seq2.length());
    
    // Forward pass
    std::vector<int> forward(n + 1, INT_MIN);
    forward[0] = 0;
    for (int j = 1; j <= n; ++j) {
        forward[j] = forward[j-1] + gap_penalty;
    }
    
    for (int i = 1; i <= mid; ++i) {
        std::vector<int> new_forward(n + 1, INT_MIN);
        new_forward[0] = i * gap_penalty;
        for (int j = 1; j <= n; ++j) {
            int s = score(seq1[i-1], seq2[j-1], match_score, mismatch_score);
            new_forward[j] = std::max({forward[j-1] + s, forward[j] + gap_penalty, new_forward[j-1] + gap_penalty});
        }
        forward = new_forward;
    }
    
    // Backward pass
    std::vector<int> backward(n + 1, INT_MIN);
    backward[n] = 0;
    for (int j = n - 1; j >= 0; --j) {
        backward[j] = backward[j+1] + gap_penalty;
    }
    
    for (int i = static_cast<int>(seq1.length()) - 1; i >= mid; --i) {
        std::vector<int> new_backward(n + 1, INT_MIN);
        new_backward[n] = (static_cast<int>(seq1.length()) - i) * gap_penalty;
        for (int j = n - 1; j >= 0; --j) {
            int s = score(seq1[i], seq2[j], match_score, mismatch_score);
            new_backward[j] = std::max({backward[j+1] + s, backward[j] + gap_penalty, new_backward[j+1] + gap_penalty});
        }
        backward = new_backward;
    }
    
    // Find best j
    int best_j = 0;
    int best_score = INT_MIN;
    for (int j = 0; j <= n; ++j) {
        int total = forward[j] + backward[j];
        if (total > best_score) {
            best_score = total;
            best_j = j;
        }
    }
    
    return best_j;
}

void AdvancedDP::reconstructAlignment(const std::string& seq1, const std::string& seq2,
                                     int mid, std::string& aligned1, std::string& aligned2) {
    // This would reconstruct the full alignment
    // Implementation depends on storing traceback information
}

