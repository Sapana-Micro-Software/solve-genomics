#include "NeedlemanWunsch.h"
#include "utils.h"
#include <algorithm>
#include <climits>
#include <vector>

AlignmentResult NeedlemanWunsch::align(const std::string& seq1, const std::string& seq2) {
    AlignmentResult result;                                                                   // Initialize alignment result
    
    if (seq1.empty() && seq2.empty()) {                                                      // Check if both sequences are empty
        return result;                                                                        // Return empty result if both empty
    }
    
    int m = seq1.length();                                                                    // Get length of first sequence
    int n = seq2.length();                                                                    // Get length of second sequence
    
    // Create scoring matrix
    std::vector<std::vector<int>> matrix(m + 1, std::vector<int>(n + 1, 0));                  // Initialize DP matrix with zeros
    
    // Initialize first row (gaps in sequence 1)
    for (int j = 1; j <= n; ++j) {                                                            // Fill first row with gap penalties
        matrix[0][j] = matrix[0][j - 1] + gap_penalty_;                                      // Cumulative gap penalty for seq1
    }
    
    // Initialize first column (gaps in sequence 2)
    for (int i = 1; i <= m; ++i) {                                                            // Fill first column with gap penalties
        matrix[i][0] = matrix[i - 1][0] + gap_penalty_;                                      // Cumulative gap penalty for seq2
    }
    
    // Fill the matrix
    for (int i = 1; i <= m; ++i) {                                                            // Iterate through first sequence
        for (int j = 1; j <= n; ++j) {                                                        // Iterate through second sequence
            // Calculate scores from three possible paths
            int match = matrix[i - 1][j - 1] + scoreMatch(seq1[i - 1], seq2[j - 1]);         // Score for match/mismatch (diagonal)
            int delete_gap = matrix[i - 1][j] + gap_penalty_;                                 // Score for gap in seq2 (up)
            int insert_gap = matrix[i][j - 1] + gap_penalty_;                                 // Score for gap in seq1 (left)
            
            // Take maximum
            matrix[i][j] = std::max({match, delete_gap, insert_gap});                         // Choose best path (global alignment)
        }
    }
    
    result.score = matrix[m][n];                                                               // Store final alignment score
    
    // Trace back from bottom-right corner
    return traceback(matrix, seq1, seq2);                                                    // Reconstruct alignment by tracing back
}

AlignmentResult NeedlemanWunsch::traceback(const std::vector<std::vector<int>>& matrix,
                                          const std::string& seq1,
                                          const std::string& seq2) {
    AlignmentResult result;                                                                   // Initialize result structure
    int m = seq1.length();                                                                    // Get length of first sequence
    int n = seq2.length();                                                                    // Get length of second sequence
    
    result.score = matrix[m][n];                                                              // Store final alignment score
    result.start_pos1 = 0;                                                                    // Global alignment starts at beginning
    result.end_pos1 = m - 1;                                                                  // Global alignment ends at end of seq1
    result.start_pos2 = 0;                                                                    // Global alignment starts at beginning
    result.end_pos2 = n - 1;                                                                  // Global alignment ends at end of seq2
    
    std::string aligned_seq1, aligned_seq2;                                                   // Strings to build aligned sequences
    int i = m, j = n;                                                                         // Start from bottom-right corner
    
    // Trace back from bottom-right to top-left
    while (i > 0 || j > 0) {                                                                  // Continue until top-left corner
        if (i > 0 && j > 0 && 
            matrix[i][j] == matrix[i - 1][j - 1] + scoreMatch(seq1[i - 1], seq2[j - 1])) {   // Check if came from diagonal
            // Match or mismatch
            aligned_seq1 = seq1[i - 1] + aligned_seq1;                                        // Add character from seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2;                                        // Add character from seq2
            i--;                                                                              // Move diagonally up-left
            j--;
        } else if (i > 0 && matrix[i][j] == matrix[i - 1][j] + gap_penalty_) {               // Check if came from above
            // Gap in sequence 2 (insertion in seq1)
            aligned_seq1 = seq1[i - 1] + aligned_seq1;                                        // Add character from seq1
            aligned_seq2 = '-' + aligned_seq2;                                               // Add gap in seq2
            i--;                                                                              // Move up
        } else if (j > 0) {                                                                   // Otherwise came from left
            // Gap in sequence 1 (insertion in seq2)
            aligned_seq1 = '-' + aligned_seq1;                                               // Add gap in seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2;                                        // Add character from seq2
            j--;                                                                              // Move left
        } else {
            // Should not happen, but break to avoid infinite loop
            break;                                                                            // Safety break to prevent infinite loop
        }
    }
    
    result.aligned_seq1 = aligned_seq1;                                                       // Store aligned first sequence
    result.aligned_seq2 = aligned_seq2;                                                       // Store aligned second sequence
    
    return result;                                                                            // Return complete alignment result
}

