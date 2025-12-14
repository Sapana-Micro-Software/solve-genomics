#ifndef ADVANCED_DP_H
#define ADVANCED_DP_H

#include "SequenceAligner.h"
#include <string>
#include <vector>

/**
 * Advanced Dynamic Programming Algorithms for Sequence Alignment
 */
class AdvancedDP {
public:
    /**
     * Space-optimized Smith-Waterman (uses only two rows)
     * Time: O(n*m), Space: O(min(n,m))
     */
    static AlignmentResult smithWatermanSpaceOptimized(const std::string& seq1,
                                                      const std::string& seq2,
                                                      int match_score = 2,
                                                      int mismatch_score = -1,
                                                      int gap_penalty = -1);
    
    /**
     * Space-optimized Needleman-Wunsch
     * Time: O(n*m), Space: O(min(n,m))
     */
    static AlignmentResult needlemanWunschSpaceOptimized(const std::string& seq1,
                                                        const std::string& seq2,
                                                        int match_score = 2,
                                                        int mismatch_score = -1,
                                                        int gap_penalty = -1);
    
    /**
     * Affine gap penalty alignment (Gotoh algorithm)
     * Uses gap opening and extension penalties
     */
    static AlignmentResult affineGapAlignment(const std::string& seq1,
                                              const std::string& seq2,
                                              int match_score = 2,
                                              int mismatch_score = -1,
                                              int gap_open = -5,
                                              int gap_extend = -1);
    
    /**
     * Banded alignment (restricts search to diagonal band)
     * Faster for similar sequences
     */
    static AlignmentResult bandedAlignment(const std::string& seq1,
                                          const std::string& seq2,
                                          int bandwidth = 10,
                                          int match_score = 2,
                                          int mismatch_score = -1,
                                          int gap_penalty = -1);
    
    /**
     * Linear space alignment (Hirschberg algorithm)
     * Time: O(n*m), Space: O(min(n,m))
     * Recursively divides problem
     */
    static AlignmentResult hirschbergAlignment(const std::string& seq1,
                                              const std::string& seq2,
                                              int match_score = 2,
                                              int mismatch_score = -1,
                                              int gap_penalty = -1);
    
private:
    /**
     * Calculate score for two characters
     */
    static int score(char a, char b, int match_score, int mismatch_score);
    
    /**
     * Hirschberg helper: find midpoint
     */
    static int findMidpoint(const std::string& seq1, const std::string& seq2,
                           int match_score, int mismatch_score, int gap_penalty);
    
    /**
     * Reconstruct alignment from midpoint
     */
    static void reconstructAlignment(const std::string& seq1, const std::string& seq2,
                                    int mid, std::string& aligned1, std::string& aligned2);
};

#endif // ADVANCED_DP_H

