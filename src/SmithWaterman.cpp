#include "SmithWaterman.h"
#include "utils.h"
#include <algorithm>
#include <climits>
#include <vector>

AlignmentResult SmithWaterman::align(const std::string& seq1, const std::string& seq2) {
    AlignmentResult result;                                                                   // Initialize alignment result
    
    if (seq1.empty() || seq2.empty()) {                                                       // Validate input sequences
        return result;                                                                        // Return empty result if invalid
    }
    
    int m = seq1.length();                                                                    // Get length of first sequence
    int n = seq2.length();                                                                    // Get length of second sequence
    
    // Create scoring matrix
    std::vector<std::vector<int>> matrix(m + 1, std::vector<int>(n + 1, 0));                  // Initialize DP matrix with zeros
    
    // Fill the matrix
    int max_score = 0;                                                                        // Track maximum alignment score
    int max_i = 0, max_j = 0;                                                                 // Track position of maximum score
    
    for (int i = 1; i <= m; ++i) {                                                            // Iterate through first sequence
        for (int j = 1; j <= n; ++j) {                                                        // Iterate through second sequence
            // Calculate scores from three possible paths
            int match = matrix[i - 1][j - 1] + scoreMatch(seq1[i - 1], seq2[j - 1]);         // Score for match/mismatch (diagonal)
            int delete_gap = matrix[i - 1][j] + gap_penalty_;                                 // Score for gap in seq2 (up)
            int insert_gap = matrix[i][j - 1] + gap_penalty_;                                 // Score for gap in seq1 (left)
            
            // Take maximum (with 0 as minimum for local alignment)
            matrix[i][j] = std::max({0, match, delete_gap, insert_gap});                     // Choose best path, minimum is 0 (local alignment)
            
            // Track maximum score position
            if (matrix[i][j] > max_score) {                                                   // Check if this is new maximum
                max_score = matrix[i][j];                                                      // Update maximum score
                max_i = i;                                                                    // Update row position
                max_j = j;                                                                    // Update column position
            }
        }
    }
    
    result.score = max_score;                                                                 // Store maximum alignment score
    
    // If no positive score found, return empty result
    if (max_score == 0) {                                                                     // Check if any alignment found
        return result;                                                                        // Return empty result if no positive score
    }
    
    // Trace back from maximum score position
    return traceback(matrix, seq1, seq2, max_i, max_j);                                      // Reconstruct alignment by tracing back
}

AlignmentResult SmithWaterman::traceback(const std::vector<std::vector<int>>& matrix,
                                        const std::string& seq1,
                                        const std::string& seq2,
                                        int max_i,
                                        int max_j) {
    AlignmentResult result;                                                                   // Initialize result structure
    result.score = matrix[max_i][max_j];                                                      // Store alignment score
    
    std::string aligned_seq1, aligned_seq2;                                                   // Strings to build aligned sequences
    int i = max_i, j = max_j;                                                                 // Start from maximum score position
    
    // Trace back until we hit a zero (local alignment boundary)
    while (i > 0 && j > 0 && matrix[i][j] > 0) {                                              // Continue until boundary or zero score
        int current = matrix[i][j];                                                            // Current cell score
        int diagonal = matrix[i - 1][j - 1];                                                  // Diagonal (match/mismatch) score
        int up = matrix[i - 1][j];                                                             // Up (gap in seq2) score
        int left = matrix[i][j - 1];                                                            // Left (gap in seq1) score
        
        if (current == diagonal + scoreMatch(seq1[i - 1], seq2[j - 1])) {                     // Check if came from diagonal
            // Match or mismatch
            aligned_seq1 = seq1[i - 1] + aligned_seq1;                                        // Add character from seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2;                                        // Add character from seq2
            i--;                                                                              // Move diagonally up-left
            j--;
        } else if (current == up + gap_penalty_) {                                            // Check if came from above
            // Gap in sequence 2 (insertion in seq1)
            aligned_seq1 = seq1[i - 1] + aligned_seq1;                                        // Add character from seq1
            aligned_seq2 = '-' + aligned_seq2;                                                // Add gap in seq2
            i--;                                                                              // Move up
        } else if (current == left + gap_penalty_) {                                           // Check if came from left
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
    result.start_pos1 = i;                                                                    // Store start position in seq1
    result.end_pos1 = max_i - 1;                                                              // Store end position in seq1
    result.start_pos2 = j;                                                                    // Store start position in seq2
    result.end_pos2 = max_j - 1;                                                              // Store end position in seq2
    
    return result;                                                                            // Return complete alignment result
}

