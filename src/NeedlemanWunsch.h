#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>

/**
 * Needleman-Wunsch Algorithm for Global Sequence Alignment
 * Aligns entire sequences from end to end
 */
class NeedlemanWunsch : public SequenceAligner {
public:
    /**
     * Align two sequences using Needleman-Wunsch algorithm
     * @param seq1 First sequence
     * @param seq2 Second sequence
     * @return AlignmentResult containing global alignment details
     */
    AlignmentResult align(const std::string& seq1, const std::string& seq2) override;
    
private:
    /**
     * Trace back from bottom-right to construct full alignment
     * @param matrix Scoring matrix
     * @param seq1 First sequence
     * @param seq2 Second sequence
     * @return AlignmentResult with aligned sequences
     */
    AlignmentResult traceback(const std::vector<std::vector<int>>& matrix,
                             const std::string& seq1,
                             const std::string& seq2);
};

#endif // NEEDLEMAN_WUNSCH_H

