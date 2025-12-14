#ifndef SEQUENCE_ALIGNER_H
#define SEQUENCE_ALIGNER_H

#include <string>
#include <vector>

/**
 * Structure to hold alignment results
 */
struct AlignmentResult {
    int score;                          // Alignment score
    std::string aligned_seq1;           // First aligned sequence (with gaps)
    std::string aligned_seq2;           // Second aligned sequence (with gaps)
    int start_pos1;                     // Start position in sequence 1
    int end_pos1;                       // End position in sequence 1
    int start_pos2;                     // Start position in sequence 2
    int end_pos2;                       // End position in sequence 2
    
    AlignmentResult() : score(0), start_pos1(0), end_pos1(0), start_pos2(0), end_pos2(0) {}
};

/**
 * Structure to hold search results (for pattern matching)
 */
struct SearchResult {
    std::vector<int> positions;          // All positions where pattern was found
    int count;                          // Number of matches found
    
    SearchResult() : count(0) {}
};

/**
 * Base class for sequence alignment algorithms
 */
class SequenceAligner {
public:
    virtual ~SequenceAligner() = default;
    
    /**
     * Align two sequences
     * @param seq1 First sequence
     * @param seq2 Second sequence
     * @return AlignmentResult containing alignment details
     */
    virtual AlignmentResult align(const std::string& seq1, const std::string& seq2) = 0;
    
    /**
     * Set scoring parameters
     * @param match_score Score for matching bases
     * @param mismatch_score Score for mismatching bases
     * @param gap_penalty Penalty for gaps (insertions/deletions)
     */
    virtual void setScoring(int match_score, int mismatch_score, int gap_penalty) {
        match_score_ = match_score;
        mismatch_score_ = mismatch_score;
        gap_penalty_ = gap_penalty;
    }
    
protected:
    int match_score_ = 2;       // Default match score
    int mismatch_score_ = -1;   // Default mismatch score
    int gap_penalty_ = -1;      // Default gap penalty
    
    /**
     * Calculate score for matching two characters
     */
    int scoreMatch(char a, char b) const {
        return (a == b) ? match_score_ : mismatch_score_;
    }
};

#endif // SEQUENCE_ALIGNER_H

