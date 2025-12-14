#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#include "SequenceAligner.h"
#include <string>
#include <vector>

/**
 * Smith-Waterman Algorithm for Local Sequence Alignment
 * Finds the best matching subsequences between two sequences
 */
class SmithWaterman : public SequenceAligner {
public:
    /**
     * Align two sequences using Smith-Waterman algorithm
     * @param seq1 First sequence
     * @param seq2 Second sequence
     * @return AlignmentResult containing local alignment details
     */
    AlignmentResult align(const std::string& seq1, const std::string& seq2) override;
    
private:
    /**
     * Trace back from maximum score position to find optimal alignment
     * @param matrix Scoring matrix
     * @param seq1 First sequence
     * @param seq2 Second sequence
     * @param max_i Row index of maximum score
     * @param max_j Column index of maximum score
     * @return AlignmentResult with aligned sequences
     */
    AlignmentResult traceback(const std::vector<std::vector<int>>& matrix,
                             const std::string& seq1,
                             const std::string& seq2,
                             int max_i,
                             int max_j);
};

#endif // SMITH_WATERMAN_H

