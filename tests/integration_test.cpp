#include <gtest/gtest.h>
#include "ExactMatch.h"
#include "NaiveSearch.h"
#include "SmithWaterman.h"
#include "NeedlemanWunsch.h"
#include "utils.h"
#include <fstream>
#include <cstdio>
#include <sstream>
#include <chrono>
#include <algorithm>

// ============================================================================
// INTEGRATION TESTS - Test components working together
// ============================================================================

class IntegrationTest : public ::testing::Test {
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

// Test complete workflow: read sequence, clean, search, align
TEST_F(IntegrationTest, CompleteWorkflow_SearchAndAlign) {
    // Simulate reading and cleaning a sequence
    std::string raw_sequence = "  atcgatcg \n\t atcg  ";
    std::string cleaned = DNAUtils::cleanSequence(raw_sequence);
    
    EXPECT_FALSE(cleaned.empty());
    EXPECT_TRUE(DNAUtils::isValidDNA(cleaned));
    
    // Search for pattern
    std::string pattern = "ATCG";
    SearchResult search_result = exact_matcher->search(cleaned, pattern);
    
    EXPECT_GT(search_result.count, 0);
    
    // Align sequences
    std::string seq2 = "ATCGATCG";
    AlignmentResult align_result = smith_waterman->align(cleaned, seq2);
    
    EXPECT_GT(align_result.score, 0);
    
    // Calculate identity
    if (!align_result.aligned_seq1.empty()) {
        double identity = DNAUtils::calculateIdentity(align_result.aligned_seq1, align_result.aligned_seq2);
        EXPECT_GE(identity, 0.0);
        EXPECT_LE(identity, 100.0);
    }
}

// Test file I/O integration
TEST_F(IntegrationTest, FileIO_ReadAndProcess) {
    // Create test file
    const char* test_file = "test_integration.txt";
    std::ofstream file(test_file);
    file << "ATCGATCGATCG";
    file.close();
    
    // Read from file
    std::string sequence = DNAUtils::readSequenceFromFile(test_file);
    EXPECT_FALSE(sequence.empty());
    EXPECT_TRUE(DNAUtils::isValidDNA(sequence));
    
    // Process the sequence
    SearchResult result = exact_matcher->search(sequence, "ATCG");
    EXPECT_GT(result.count, 0);
    
    // Clean up
    std::remove(test_file);
}

// Test multiple algorithms on same data
TEST_F(IntegrationTest, MultipleAlgorithms_SameData) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    // Use all algorithms
    SearchResult exact_result = exact_matcher->search(seq1, "ATCG");
    SearchResult naive_result = naive_searcher->search(seq1, "ATCG");
    AlignmentResult sw_result = smith_waterman->align(seq1, seq2);
    AlignmentResult nw_result = needleman_wunsch->align(seq1, seq2);
    
    // All should produce valid results
    EXPECT_EQ(exact_result.count, naive_result.count);
    EXPECT_GT(sw_result.score, 0);
    EXPECT_GT(nw_result.score, 0);
}

// Test error handling across components
TEST_F(IntegrationTest, ErrorHandling_InvalidInputs) {
    // Invalid DNA sequence
    std::string invalid = "ATCGX";
    std::string cleaned = DNAUtils::cleanSequence(invalid);
    EXPECT_TRUE(cleaned.empty());
    
    // Empty sequences
    SearchResult search_result = exact_matcher->search("", "ATCG");
    EXPECT_EQ(search_result.count, 0);
    
    AlignmentResult align_result = smith_waterman->align("", "ATCG");
    EXPECT_EQ(align_result.score, 0);
}

// ============================================================================
// UX TESTS - Test user experience and usability
// ============================================================================

class UXTest : public ::testing::Test {
protected:
    void SetUp() override {
        exact_matcher = new ExactMatch();
        smith_waterman = new SmithWaterman();
        smith_waterman->setScoring(2, -1, -1);
    }
    
    void TearDown() override {
        delete exact_matcher;
        delete smith_waterman;
    }
    
    ExactMatch* exact_matcher;
    SmithWaterman* smith_waterman;
};

// Test that results are consistent and predictable
TEST_F(UXTest, Consistency_SameInputSameOutput) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    // Run multiple times
    AlignmentResult result1 = smith_waterman->align(seq1, seq2);
    AlignmentResult result2 = smith_waterman->align(seq1, seq2);
    AlignmentResult result3 = smith_waterman->align(seq1, seq2);
    
    // Should get same results
    EXPECT_EQ(result1.score, result2.score);
    EXPECT_EQ(result2.score, result3.score);
}

// Test that algorithms handle reasonable input sizes
TEST_F(UXTest, Performance_ReasonableInputSizes) {
    // Create sequences of reasonable size
    std::string seq1(100, 'A');
    std::string seq2(100, 'A');
    
    auto start = std::chrono::high_resolution_clock::now();
    AlignmentResult result = smith_waterman->align(seq1, seq2);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Should complete in reasonable time (< 1 second)
    EXPECT_LT(duration.count(), 1000);
    EXPECT_GT(result.score, 0);
}

// Test that error messages/behavior are clear
TEST_F(UXTest, ErrorHandling_ClearBehavior) {
    // Empty input should return empty/zero results, not crash
    SearchResult search_result = exact_matcher->search("", "ATCG");
    EXPECT_EQ(search_result.count, 0);
    EXPECT_TRUE(search_result.positions.empty());
    
    AlignmentResult align_result = smith_waterman->align("", "ATCG");
    EXPECT_EQ(align_result.score, 0);
    EXPECT_TRUE(align_result.aligned_seq1.empty());
}

// Test that results are interpretable
TEST_F(UXTest, Interpretability_ResultsMakeSense) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCG";
    
    AlignmentResult result = smith_waterman->align(seq1, seq2);
    
    // Results should be interpretable
    EXPECT_GT(result.score, 0);  // Positive score for match
    EXPECT_FALSE(result.aligned_seq1.empty());  // Should have alignment
    EXPECT_EQ(result.aligned_seq1.length(), result.aligned_seq2.length());  // Same length
    
    // Identity should be calculable
    if (!result.aligned_seq1.empty()) {
        double identity = DNAUtils::calculateIdentity(result.aligned_seq1, result.aligned_seq2);
        EXPECT_GE(identity, 0.0);
        EXPECT_LE(identity, 100.0);
    }
}

// Test case-insensitive behavior (user-friendly)
TEST_F(UXTest, UserFriendly_CaseInsensitive) {
    std::string seq1 = "atcg";
    std::string seq2 = "ATCG";
    
    SearchResult result = exact_matcher->search(seq1, seq2);
    EXPECT_GT(result.count, 0);  // Should find matches despite case difference
}

// Test that formatting works correctly
TEST_F(UXTest, Formatting_AlignmentDisplay) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    
    AlignmentResult result = smith_waterman->align(seq1, seq2);
    
    if (!result.aligned_seq1.empty()) {
        std::string formatted = DNAUtils::formatAlignment(result.aligned_seq1, result.aligned_seq2);
        
        // Should contain expected elements
        EXPECT_NE(formatted.find("Seq1:"), std::string::npos);
        EXPECT_NE(formatted.find("Seq2:"), std::string::npos);
    }
}

// Test edge cases that users might encounter
TEST_F(UXTest, EdgeCases_CommonUserScenarios) {
    // Single character
    SearchResult result1 = exact_matcher->search("A", "A");
    EXPECT_EQ(result1.count, 1);
    
    // Pattern longer than sequence
    SearchResult result2 = exact_matcher->search("ATC", "ATCG");
    EXPECT_EQ(result2.count, 0);
    
    // Very long pattern
    std::string long_seq(1000, 'A');
    std::string pattern(100, 'A');
    SearchResult result3 = exact_matcher->search(long_seq, pattern);
    EXPECT_GT(result3.count, 0);
}

// Test that scoring parameters are configurable
TEST_F(UXTest, Configurability_ScoringParameters) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCG";
    
    // Default scoring
    AlignmentResult result1 = smith_waterman->align(seq1, seq2);
    
    // Custom scoring
    smith_waterman->setScoring(10, -5, -3);
    AlignmentResult result2 = smith_waterman->align(seq1, seq2);
    
    // Should produce different scores
    EXPECT_NE(result1.score, result2.score);
    
    // Reset
    smith_waterman->setScoring(2, -1, -1);
}

// Test that algorithms work with real-world-like data
TEST_F(UXTest, RealWorld_RealisticSequences) {
    // Simulate a gene sequence
    std::string gene1 = "ATCGATCGATCGATCGATCG";
    std::string gene2 = "ATCGATCGATCGATCGATCC";  // One mutation
    
    AlignmentResult result = smith_waterman->align(gene1, gene2);
    
    // Should find good alignment despite mutation
    EXPECT_GT(result.score, 0);
    
    // Calculate similarity
    if (!result.aligned_seq1.empty()) {
        double identity = DNAUtils::calculateIdentity(result.aligned_seq1, result.aligned_seq2);
        EXPECT_GT(identity, 80.0);  // Should be highly similar
    }
}

// Test memory safety (no crashes on large inputs)
TEST_F(UXTest, Stability_LargeInputs) {
    // Create moderately large sequences
    std::string seq1(500, 'A');
    std::string seq2(500, 'A');
    
    // Should not crash
    EXPECT_NO_THROW({
        AlignmentResult result = smith_waterman->align(seq1, seq2);
        EXPECT_GT(result.score, 0);
    });
}

// Test that all algorithms can be used together
TEST_F(UXTest, Compatibility_AllAlgorithmsTogether) {
    std::string seq = "ATCGATCGATCG";
    std::string pattern = "ATCG";
    std::string seq2 = "ATCGATCGATCG";
    
    NaiveSearch* naive_searcher = new NaiveSearch();
    
    // Use all algorithms
    SearchResult exact = exact_matcher->search(seq, pattern);
    SearchResult naive = naive_searcher->search(seq, pattern);
    AlignmentResult sw = smith_waterman->align(seq, seq2);
    NeedlemanWunsch nw_aligner;
    AlignmentResult nw = nw_aligner.align(seq, seq2);
    
    // All should work
    EXPECT_GT(exact.count, 0);
    EXPECT_GT(naive.count, 0);
    EXPECT_GT(sw.score, 0);
    EXPECT_GT(nw.score, 0);
    
    delete naive_searcher;
}

