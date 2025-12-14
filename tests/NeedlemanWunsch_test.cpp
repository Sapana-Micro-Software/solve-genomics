#include <gtest/gtest.h>
#include "NeedlemanWunsch.h"
#include "utils.h"
#include <algorithm>

// ============================================================================
// UNIT TESTS - Testing individual methods
// ============================================================================

class NeedlemanWunschUnitTest : public ::testing::Test {
protected:
    void SetUp() override {
        aligner = new NeedlemanWunsch();
        aligner->setScoring(2, -1, -1);  // match=2, mismatch=-1, gap=-1
    }
    
    void TearDown() override {
        delete aligner;
    }
    
    NeedlemanWunsch* aligner;
};

TEST_F(NeedlemanWunschUnitTest, Align_IdenticalSequences) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCG";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    EXPECT_GT(result.score, 0);
    EXPECT_FALSE(result.aligned_seq1.empty());
    EXPECT_FALSE(result.aligned_seq2.empty());
    EXPECT_EQ(result.aligned_seq1.length(), result.aligned_seq2.length());
    // Perfect match should have high score
    EXPECT_GE(result.score, 8);  // 4 matches * 2 = 8
}

TEST_F(NeedlemanWunschUnitTest, Align_SequencesWithMismatch) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCC";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    EXPECT_GT(result.score, 0);
    // 3 matches * 2 + 1 mismatch * (-1) = 6 - 1 = 5
    EXPECT_GE(result.score, 5);
}

TEST_F(NeedlemanWunschUnitTest, Align_SequencesWithGaps) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATC";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    EXPECT_GT(result.score, 0);
    EXPECT_FALSE(result.aligned_seq1.empty());
    EXPECT_FALSE(result.aligned_seq2.empty());
    // Should have a gap in one sequence
    EXPECT_TRUE(result.aligned_seq1.find('-') != std::string::npos || 
                result.aligned_seq2.find('-') != std::string::npos);
}

TEST_F(NeedlemanWunschUnitTest, Align_EmptyInput) {
    AlignmentResult result1 = aligner->align("", "ATCG");
    // Empty vs non-empty should have negative score (all gaps)
    EXPECT_LE(result1.score, 0);
    
    AlignmentResult result2 = aligner->align("ATCG", "");
    EXPECT_LE(result2.score, 0);
    
    AlignmentResult result3 = aligner->align("", "");
    EXPECT_EQ(result3.score, 0);
}

TEST_F(NeedlemanWunschUnitTest, Align_DifferentLengthSequences) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    EXPECT_FALSE(result.aligned_seq1.empty());
    EXPECT_FALSE(result.aligned_seq2.empty());
    EXPECT_EQ(result.aligned_seq1.length(), result.aligned_seq2.length());
    // Should align entire sequences with gaps
}

TEST_F(NeedlemanWunschUnitTest, SetScoring_CustomScores) {
    aligner->setScoring(5, -2, -3);
    
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCG";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    // With higher match score, should get higher total score
    EXPECT_GT(result.score, 8);  // Should be higher than default scoring
}

TEST_F(NeedlemanWunschUnitTest, Align_GlobalAlignmentProperty) {
    // Global alignment should align entire sequences
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCG";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    EXPECT_EQ(result.start_pos1, 0);
    EXPECT_EQ(result.end_pos1, static_cast<int>(seq1.length()) - 1);
    EXPECT_EQ(result.start_pos2, 0);
    EXPECT_EQ(result.end_pos2, static_cast<int>(seq2.length()) - 1);
}

// ============================================================================
// REGRESSION TESTS - Test known good outputs
// ============================================================================

TEST_F(NeedlemanWunschUnitTest, Regression_KnownAlignment1) {
    // Known test case: identical sequences
    std::string seq1 = "ACGTACGT";
    std::string seq2 = "ACGTACGT";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    // Perfect match: 8 matches * 2 = 16
    EXPECT_EQ(result.score, 16);
    EXPECT_EQ(result.aligned_seq1, seq1);
    EXPECT_EQ(result.aligned_seq2, seq2);
}

TEST_F(NeedlemanWunschUnitTest, Regression_KnownAlignment2) {
    // Known test case: one gap
    std::string seq1 = "ATCG";
    std::string seq2 = "ATC";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    // 3 matches * 2 + 1 gap * (-1) = 6 - 1 = 5
    EXPECT_EQ(result.score, 5);
    EXPECT_EQ(result.aligned_seq1.length(), result.aligned_seq2.length());
}

TEST_F(NeedlemanWunschUnitTest, Regression_ScoringConsistency) {
    // Same sequences should give same score
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult result1 = aligner->align(seq1, seq2);
    AlignmentResult result2 = aligner->align(seq1, seq2);
    
    EXPECT_EQ(result1.score, result2.score);
}

// ============================================================================
// BLACKBOX TESTS - Test without knowing implementation
// ============================================================================

TEST_F(NeedlemanWunschUnitTest, Blackbox_AlignedSequencesSameLength) {
    // Aligned sequences should always have the same length
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATCG", "ATCG"},
        {"ATCGATCG", "ATCG"},
        {"ATCG", "ATCGATCG"},
        {"AAAA", "TTTT"},
        {"A", "ATCG"}
    };
    
    for (const auto& test : test_cases) {
        AlignmentResult result = aligner->align(test.first, test.second);
        if (!result.aligned_seq1.empty()) {
            EXPECT_EQ(result.aligned_seq1.length(), result.aligned_seq2.length());
        }
    }
}

TEST_F(NeedlemanWunschUnitTest, Blackbox_GlobalAlignmentProperty) {
    // Global alignment should always align from start to end
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    EXPECT_EQ(result.start_pos1, 0);
    EXPECT_EQ(result.start_pos2, 0);
    EXPECT_EQ(result.end_pos1, static_cast<int>(seq1.length()) - 1);
    EXPECT_EQ(result.end_pos2, static_cast<int>(seq2.length()) - 1);
}

TEST_F(NeedlemanWunschUnitTest, Blackbox_CommutativeProperty) {
    // Aligning seq1 vs seq2 should give same score as seq2 vs seq1
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult result1 = aligner->align(seq1, seq2);
    AlignmentResult result2 = aligner->align(seq2, seq1);
    
    EXPECT_EQ(result1.score, result2.score);
}

TEST_F(NeedlemanWunschUnitTest, Blackbox_IdentityCalculation) {
    // Calculate identity should work on aligned sequences
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    if (!result.aligned_seq1.empty()) {
        double identity = DNAUtils::calculateIdentity(result.aligned_seq1, result.aligned_seq2);
        EXPECT_GE(identity, 0.0);
        EXPECT_LE(identity, 100.0);
        
        // For identical sequences, identity should be 100%
        EXPECT_DOUBLE_EQ(identity, 100.0);
    }
}

TEST_F(NeedlemanWunschUnitTest, Blackbox_ScoreIncreasesWithSimilarity) {
    // More similar sequences should have higher scores
    std::string base = "ATCGATCG";
    
    std::string identical = "ATCGATCG";
    std::string similar = "ATCGATCC";
    std::string different = "GGGGGGGG";
    
    AlignmentResult result1 = aligner->align(base, identical);
    AlignmentResult result2 = aligner->align(base, similar);
    AlignmentResult result3 = aligner->align(base, different);
    
    EXPECT_GE(result1.score, result2.score);
    EXPECT_GE(result2.score, result3.score);
}

TEST_F(NeedlemanWunschUnitTest, Blackbox_AlignedSequenceContainsOriginal) {
    // The aligned sequence should contain all characters from original (or gaps)
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCG";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    if (!result.aligned_seq1.empty()) {
        // Remove gaps
        std::string seq1_no_gaps = result.aligned_seq1;
        seq1_no_gaps.erase(std::remove(seq1_no_gaps.begin(), seq1_no_gaps.end(), '-'), 
                          seq1_no_gaps.end());
        
        // Should contain all characters from original (in order for global alignment)
        EXPECT_EQ(seq1_no_gaps.length(), seq1.length());
    }
}

TEST_F(NeedlemanWunschUnitTest, Blackbox_EmptySequenceHandling) {
    // Test various empty sequence combinations
    AlignmentResult result1 = aligner->align("", "");
    EXPECT_EQ(result1.score, 0);
    
    AlignmentResult result2 = aligner->align("ATCG", "");
    EXPECT_LE(result2.score, 0);  // All gaps = negative score
    
    AlignmentResult result3 = aligner->align("", "ATCG");
    EXPECT_LE(result3.score, 0);  // All gaps = negative score
}

TEST_F(NeedlemanWunschUnitTest, Blackbox_LengthConservation) {
    // Global alignment: sum of non-gap characters should equal original length
    std::string seq1 = "ATCG";
    std::string seq2 = "ATC";
    
    AlignmentResult result = aligner->align(seq1, seq2);
    
    if (!result.aligned_seq1.empty()) {
        std::string seq1_no_gaps = result.aligned_seq1;
        seq1_no_gaps.erase(std::remove(seq1_no_gaps.begin(), seq1_no_gaps.end(), '-'), 
                          seq1_no_gaps.end());
        EXPECT_EQ(seq1_no_gaps.length(), seq1.length());
        
        std::string seq2_no_gaps = result.aligned_seq2;
        seq2_no_gaps.erase(std::remove(seq2_no_gaps.begin(), seq2_no_gaps.end(), '-'), 
                          seq2_no_gaps.end());
        EXPECT_EQ(seq2_no_gaps.length(), seq2.length());
    }
}

