#include <gtest/gtest.h>
#include "ExactMatch.h"
#include "NaiveSearch.h"
#include "SmithWaterman.h"
#include "NeedlemanWunsch.h"
#include "utils.h"
#include <chrono>
#include <algorithm>

// ============================================================================
// A-B COMPARISON TESTS - Compare different algorithms on same inputs
// ============================================================================

class AlgorithmComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        exact_matcher = new ExactMatch();
        naive_searcher = new NaiveSearch();
        smith_waterman = new SmithWaterman();
        needleman_wunsch = new NeedlemanWunsch();
        
        smith_waterman->setScoring(2, -1, -1);
        needleman_wunsch->setScoring(2, -1, -1);
    }
    
    void TearDown() override {
        delete exact_matcher;
        delete naive_searcher;
        delete smith_waterman;
        delete needleman_wunsch;
    }
    
    ExactMatch* exact_matcher;
    NaiveSearch* naive_searcher;
    SmithWaterman* smith_waterman;
    NeedlemanWunsch* needleman_wunsch;
};

// Compare ExactMatch vs NaiveSearch (should produce identical results)
TEST_F(AlgorithmComparisonTest, AB_ExactMatchVsNaiveSearch) {
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATCGATCG", "ATCG"},
        {"AAAAA", "AAA"},
        {"ATCGATCGATCG", "TCGA"},
        {"ATCG", "GGGG"},
        {"A", "A"}
    };
    
    for (const auto& test : test_cases) {
        SearchResult exact_result = exact_matcher->search(test.first, test.second);
        SearchResult naive_result = naive_searcher->search(test.first, test.second);
        
        EXPECT_EQ(exact_result.count, naive_result.count)
            << "Different counts for sequence: " << test.first 
            << ", pattern: " << test.second;
        
        EXPECT_EQ(exact_result.positions.size(), naive_result.positions.size());
        
        for (size_t i = 0; i < exact_result.positions.size(); ++i) {
            EXPECT_EQ(exact_result.positions[i], naive_result.positions[i])
                << "Different position at index " << i;
        }
    }
}

// Compare Smith-Waterman vs Needleman-Wunsch on identical sequences
TEST_F(AlgorithmComparisonTest, AB_SmithWatermanVsNeedlemanWunsch_Identical) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    
    // For identical sequences, both should find perfect alignment
    // Needleman-Wunsch should have higher or equal score (global alignment)
    EXPECT_GE(nw_result.score, sw_result.score);
    
    // Both should have high scores
    EXPECT_GT(sw_result.score, 0);
    EXPECT_GT(nw_result.score, 0);
}

// Compare Smith-Waterman vs Needleman-Wunsch on similar sequences
TEST_F(AlgorithmComparisonTest, AB_SmithWatermanVsNeedlemanWunsch_Similar) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCC";
    
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    
    // Both should find good alignments
    EXPECT_GT(sw_result.score, 0);
    EXPECT_GT(nw_result.score, 0);
    
    // Needleman-Wunsch aligns globally, so it might have different score
    // Smith-Waterman finds local best, which might be better for local similarity
}

// Compare Smith-Waterman vs Needleman-Wunsch on sequences with only local similarity
TEST_F(AlgorithmComparisonTest, AB_SmithWatermanVsNeedlemanWunsch_LocalSimilarity) {
    std::string seq1 = "XXXXXATCGATCGYYYYY";
    std::string seq2 = "ZZZZZATCGATCGWWWWW";
    
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    
    // Smith-Waterman should find the local similarity
    EXPECT_GT(sw_result.score, 0);
    
    // Needleman-Wunsch will align globally, including the different regions
    EXPECT_GT(nw_result.score, 0);
    
    // Smith-Waterman should have better score for local similarity
    EXPECT_GE(sw_result.score, nw_result.score);
}

// Compare alignment algorithms on different sequence lengths
TEST_F(AlgorithmComparisonTest, AB_DifferentSequenceLengths) {
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATCG", "ATCGATCG"},      // seq2 longer
        {"ATCGATCG", "ATCG"},      // seq1 longer
        {"ATCG", "ATCG"},          // same length
        {"A", "ATCGATCG"},         // very different lengths
        {"ATCGATCG", "A"}          // very different lengths
    };
    
    for (const auto& test : test_cases) {
        AlignmentResult sw_result = smith_waterman->align(test.first, test.second);
        AlignmentResult nw_result = needleman_wunsch->align(test.first, test.second);
        
        // Both should produce valid results
        EXPECT_GE(sw_result.score, 0);
        EXPECT_GE(nw_result.score, -100);  // Can be negative for global with many gaps
        
        // Aligned sequences should have same length
        if (!sw_result.aligned_seq1.empty()) {
            EXPECT_EQ(sw_result.aligned_seq1.length(), sw_result.aligned_seq2.length());
        }
        if (!nw_result.aligned_seq1.empty()) {
            EXPECT_EQ(nw_result.aligned_seq1.length(), nw_result.aligned_seq2.length());
        }
    }
}

// Performance comparison (A-B testing for performance)
TEST_F(AlgorithmComparisonTest, AB_PerformanceComparison) {
    std::string seq1 = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    std::string seq2 = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    
    // Time Smith-Waterman
    auto start_sw = std::chrono::high_resolution_clock::now();
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    auto end_sw = std::chrono::high_resolution_clock::now();
    auto duration_sw = std::chrono::duration_cast<std::chrono::microseconds>(end_sw - start_sw);
    
    // Time Needleman-Wunsch
    auto start_nw = std::chrono::high_resolution_clock::now();
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    auto end_nw = std::chrono::high_resolution_clock::now();
    auto duration_nw = std::chrono::duration_cast<std::chrono::microseconds>(end_nw - start_nw);
    
    // Both should complete (no timeout)
    EXPECT_GT(sw_result.score, 0);
    EXPECT_GT(nw_result.score, 0);
    
    // Both should complete in reasonable time (< 1 second for small sequences)
    EXPECT_LT(duration_sw.count(), 1000000);
    EXPECT_LT(duration_nw.count(), 1000000);
    
    // Performance should be similar (both O(n*m))
    // Allow 10x difference for implementation variations
    EXPECT_LT(std::max(duration_sw.count(), duration_nw.count()) / 
              std::min(duration_sw.count(), duration_nw.count()), 10);
}

// Compare scoring systems
TEST_F(AlgorithmComparisonTest, AB_DifferentScoringSystems) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCG";
    
    // Default scoring
    AlignmentResult result1 = smith_waterman->align(seq1, seq2);
    
    // Higher match score
    smith_waterman->setScoring(5, -1, -1);
    AlignmentResult result2 = smith_waterman->align(seq1, seq2);
    
    // Higher match score should give higher total score
    EXPECT_GT(result2.score, result1.score);
    
    // Reset
    smith_waterman->setScoring(2, -1, -1);
}

// Compare on edge cases
TEST_F(AlgorithmComparisonTest, AB_EdgeCases) {
    // Single character
    AlignmentResult sw1 = smith_waterman->align("A", "A");
    AlignmentResult nw1 = needleman_wunsch->align("A", "A");
    EXPECT_GT(sw1.score, 0);
    EXPECT_GT(nw1.score, 0);
    
    // Single character mismatch
    AlignmentResult sw2 = smith_waterman->align("A", "T");
    AlignmentResult nw2 = needleman_wunsch->align("A", "T");
    EXPECT_GE(sw2.score, 0);  // Local might be 0
    EXPECT_LE(nw2.score, 0);  // Global will have mismatch penalty
    
    // One empty
    AlignmentResult sw3 = smith_waterman->align("ATCG", "");
    AlignmentResult nw3 = needleman_wunsch->align("ATCG", "");
    EXPECT_EQ(sw3.score, 0);  // Local: no positive alignment
    EXPECT_LE(nw3.score, 0);  // Global: all gaps
}

// Compare identity calculations
TEST_F(AlgorithmComparisonTest, AB_IdentityCalculations) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    
    if (!sw_result.aligned_seq1.empty() && !nw_result.aligned_seq1.empty()) {
        double sw_identity = DNAUtils::calculateIdentity(sw_result.aligned_seq1, sw_result.aligned_seq2);
        double nw_identity = DNAUtils::calculateIdentity(nw_result.aligned_seq1, nw_result.aligned_seq2);
        
        // For identical sequences, both should have high identity
        EXPECT_GT(sw_identity, 90.0);
        EXPECT_GT(nw_identity, 90.0);
    }
}

// Compare on sequences with gaps
TEST_F(AlgorithmComparisonTest, AB_SequencesWithGaps) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATC";
    
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    
    // Both should handle gaps
    EXPECT_GT(sw_result.score, 0);
    EXPECT_GT(nw_result.score, 0);
    
    // Both should produce aligned sequences with same length
    if (!sw_result.aligned_seq1.empty()) {
        EXPECT_EQ(sw_result.aligned_seq1.length(), sw_result.aligned_seq2.length());
    }
    if (!nw_result.aligned_seq1.empty()) {
        EXPECT_EQ(nw_result.aligned_seq1.length(), nw_result.aligned_seq2.length());
    }
}

// Regression: Compare known test cases
TEST_F(AlgorithmComparisonTest, AB_Regression_KnownTestCases) {
    // Known test case from literature
    std::string seq1 = "ACGTACGT";
    std::string seq2 = "ACGTACGT";
    
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    
    // Both should find perfect alignment
    EXPECT_GT(sw_result.score, 10);
    EXPECT_GT(nw_result.score, 10);
    
    // Needleman-Wunsch should align globally (full sequence)
    EXPECT_EQ(nw_result.start_pos1, 0);
    EXPECT_EQ(nw_result.end_pos1, static_cast<int>(seq1.length()) - 1);
}

