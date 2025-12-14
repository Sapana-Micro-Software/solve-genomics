#ifndef WARP_CTC_H
#define WARP_CTC_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>

/**
 * WARP-CTC (Connectionist Temporal Classification) for DNA pattern matching
 * Handles sequence alignment with insertions, deletions, and substitutions
 * Useful for finding patterns in sequences with gaps and variations
 */
class WarpCTC {
public:
    struct CTCResult {
        std::vector<size_t> positions;          // Positions where pattern aligns
        std::vector<double> alignment_scores;   // Alignment scores at each position
        std::string aligned_pattern;            // Aligned pattern (with gaps)
        std::string aligned_sequence;           // Aligned sequence region
        double total_score;                     // Total CTC score
        int num_matches;                        // Number of matches found
        
        CTCResult() : total_score(0.0), num_matches(0) {}
    };
    
    WarpCTC(double blank_probability = 0.1);
    
    /**
     * Find pattern in sequence using CTC alignment
     * @param sequence DNA sequence to search in
     * @param pattern Pattern to find
     * @param max_gaps Maximum number of gaps allowed
     * @return CTCResult with alignments
     */
    CTCResult search(const std::string& sequence, 
                    const std::string& pattern,
                    int max_gaps = 5);
    
    /**
     * Calculate CTC loss/score for alignment
     * @param sequence Sequence
     * @param pattern Pattern
     * @return CTC score (higher is better)
     */
    double calculateCTCLoss(const std::string& sequence, const std::string& pattern);
    
    /**
     * Forward algorithm for CTC
     * Computes probability of sequence given pattern
     */
    double forwardAlgorithm(const std::string& sequence, const std::string& pattern);
    
    /**
     * Backward algorithm for CTC
     * Computes backward probabilities
     */
    double backwardAlgorithm(const std::string& sequence, const std::string& pattern);
    
    /**
     * Viterbi decoding (best path)
     * Finds most likely alignment path
     */
    std::string viterbiDecode(const std::string& sequence, const std::string& pattern);
    
    /**
     * Beam search for pattern matching
     * @param sequence Sequence to search
     * @param pattern Pattern to find
     * @param beam_width Beam width
     * @return Best alignment paths
     */
    std::vector<std::pair<std::string, double>> beamSearch(const std::string& sequence,
                                                           const std::string& pattern,
                                                           int beam_width = 5);
    
    /**
     * Align pattern to sequence with CTC
     * @param sequence Sequence
     * @param pattern Pattern
     * @return Aligned pair (pattern, sequence)
     */
    std::pair<std::string, std::string> align(const std::string& sequence,
                                              const std::string& pattern);
    
private:
    double blank_probability_;
    
    /**
     * Create extended pattern with blanks (CTC format)
     * Pattern "ATCG" becomes " A T C G " (with blanks)
     */
    std::string extendWithBlanks(const std::string& pattern);
    
    /**
     * Calculate emission probability
     * Probability of observing sequence character given pattern character
     */
    double emissionProbability(char seq_char, char pattern_char);
    
    /**
     * Calculate transition probability in CTC lattice
     */
    double transitionProbability(int from_state, int to_state, const std::string& extended_pattern);
    
    /**
     * Build CTC alignment matrix
     */
    std::vector<std::vector<double>> buildAlignmentMatrix(const std::string& sequence,
                                                          const std::string& extended_pattern);
    
    /**
     * Trace back to find best alignment
     */
    std::pair<std::string, std::string> tracebackAlignment(
        const std::vector<std::vector<double>>& matrix,
        const std::string& sequence,
        const std::string& extended_pattern);
    
    /**
     * Find all valid alignments using CTC
     */
    std::vector<std::pair<size_t, double>> findAllAlignments(
        const std::string& sequence,
        const std::string& pattern);
};

#endif // WARP_CTC_H

