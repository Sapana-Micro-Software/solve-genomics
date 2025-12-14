#include <gtest/gtest.h>
#include "utils.h"
#include <fstream>
#include <cstdio>
#include <algorithm>

using namespace DNAUtils;

// ============================================================================
// UNIT TESTS - Testing individual utility functions
// ============================================================================

class UtilsUnitTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilsUnitTest, IsValidDNA_ValidSequence) {
    EXPECT_TRUE(isValidDNA("ATCG"));
    EXPECT_TRUE(isValidDNA("atcg"));
    EXPECT_TRUE(isValidDNA("ATCGATCG"));
    EXPECT_TRUE(isValidDNA("A"));
    EXPECT_TRUE(isValidDNA("T"));
    EXPECT_TRUE(isValidDNA("G"));
    EXPECT_TRUE(isValidDNA("C"));
}

TEST_F(UtilsUnitTest, IsValidDNA_InvalidSequence) {
    EXPECT_FALSE(isValidDNA(""));
    EXPECT_FALSE(isValidDNA("ATCGX"));
    EXPECT_FALSE(isValidDNA("123"));
    EXPECT_FALSE(isValidDNA("ATCG "));
    EXPECT_FALSE(isValidDNA("ATCG\n"));
}

TEST_F(UtilsUnitTest, IsValidDNA_WithGaps) {
    EXPECT_TRUE(isValidDNA("AT-CG"));
    EXPECT_TRUE(isValidDNA("A-T-C-G"));
    EXPECT_TRUE(isValidDNA("---"));
}

TEST_F(UtilsUnitTest, ToUpperCase) {
    EXPECT_EQ(toUpperCase("atcg"), "ATCG");
    EXPECT_EQ(toUpperCase("ATCG"), "ATCG");
    EXPECT_EQ(toUpperCase("AtCg"), "ATCG");
    EXPECT_EQ(toUpperCase(""), "");
}

TEST_F(UtilsUnitTest, RemoveWhitespace) {
    EXPECT_EQ(removeWhitespace("ATCG"), "ATCG");
    EXPECT_EQ(removeWhitespace("A T C G"), "ATCG");
    EXPECT_EQ(removeWhitespace("A\tT\nC\rG"), "ATCG");
    EXPECT_EQ(removeWhitespace("  ATCG  "), "ATCG");
    EXPECT_EQ(removeWhitespace(""), "");
}

TEST_F(UtilsUnitTest, CleanSequence_ValidInput) {
    std::string cleaned = cleanSequence("atcg");
    EXPECT_EQ(cleaned, "ATCG");
    
    cleaned = cleanSequence("A T C G");
    EXPECT_EQ(cleaned, "ATCG");
    
    cleaned = cleanSequence("  atcg  ");
    EXPECT_EQ(cleaned, "ATCG");
}

TEST_F(UtilsUnitTest, CleanSequence_InvalidInput) {
    EXPECT_EQ(cleanSequence("ATCGX"), "");
    EXPECT_EQ(cleanSequence("123"), "");
    EXPECT_EQ(cleanSequence(""), "");
}

TEST_F(UtilsUnitTest, ReadSequenceFromFile_ValidFile) {
    // Create a temporary test file
    const char* test_file = "test_sequence.txt";
    std::ofstream file(test_file);
    file << "ATCGATCG";
    file.close();
    
    std::string sequence = readSequenceFromFile(test_file);
    EXPECT_EQ(sequence, "ATCGATCG");
    
    // Clean up
    std::remove(test_file);
}

TEST_F(UtilsUnitTest, ReadSequenceFromFile_InvalidFile) {
    std::string sequence = readSequenceFromFile("nonexistent_file.txt");
    EXPECT_EQ(sequence, "");
}

TEST_F(UtilsUnitTest, FormatAlignment_ValidAlignment) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    std::string formatted = formatAlignment(seq1, seq2, 4);
    
    EXPECT_NE(formatted.find("Seq1:"), std::string::npos);
    EXPECT_NE(formatted.find("Seq2:"), std::string::npos);
    EXPECT_NE(formatted.find("ATCG"), std::string::npos);
}

TEST_F(UtilsUnitTest, FormatAlignment_MismatchedLength) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCGATCG";
    std::string formatted = formatAlignment(seq1, seq2);
    
    EXPECT_NE(formatted.find("Error"), std::string::npos);
}

TEST_F(UtilsUnitTest, CalculateIdentity_PerfectMatch) {
    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    double identity = calculateIdentity(seq1, seq2);
    EXPECT_DOUBLE_EQ(identity, 100.0);
}

TEST_F(UtilsUnitTest, CalculateIdentity_PartialMatch) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCC";
    double identity = calculateIdentity(seq1, seq2);
    EXPECT_DOUBLE_EQ(identity, 75.0);  // 3 out of 4 match
}

TEST_F(UtilsUnitTest, CalculateIdentity_WithGaps) {
    std::string seq1 = "AT-CG";
    std::string seq2 = "ATCCG";
    double identity = calculateIdentity(seq1, seq2);
    // Should calculate based on non-gap positions
    EXPECT_GT(identity, 0.0);
    EXPECT_LE(identity, 100.0);
}

TEST_F(UtilsUnitTest, CalculateIdentity_MismatchedLength) {
    std::string seq1 = "ATCG";
    std::string seq2 = "ATCGATCG";
    double identity = calculateIdentity(seq1, seq2);
    EXPECT_DOUBLE_EQ(identity, 0.0);
}

// ============================================================================
// REGRESSION TESTS - Test known good outputs
// ============================================================================

TEST_F(UtilsUnitTest, Regression_CleanSequence_ComplexInput) {
    // Known good output for complex input
    std::string input = "  atcg \n\t atcg  ";
    std::string expected = "ATCGATCG";
    std::string result = cleanSequence(input);
    EXPECT_EQ(result, expected);
}

// ============================================================================
// BLACKBOX TESTS - Test without knowing implementation
// ============================================================================

TEST_F(UtilsUnitTest, Blackbox_CleanSequence_VariousInputs) {
    // Test various inputs and verify output is always valid DNA or empty
    std::vector<std::string> inputs = {
        "ATCG", "atcg", "A T C G", "  ATCG  ", "ATCGX", "123", ""
    };
    
    for (const auto& input : inputs) {
        std::string result = cleanSequence(input);
        if (!result.empty()) {
            // If not empty, should be valid DNA
            std::string for_validation = result;
            for_validation.erase(
                std::remove(for_validation.begin(), for_validation.end(), '-'),
                for_validation.end()
            );
            EXPECT_TRUE(isValidDNA(for_validation));
        }
    }
}

TEST_F(UtilsUnitTest, Blackbox_CalculateIdentity_Properties) {
    // Identity should always be between 0 and 100
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATCG", "ATCG"},
        {"ATCG", "AAAA"},
        {"ATCG", "TTTT"},
        {"AT-CG", "ATCCG"},
        {"", ""}
    };
    
    for (const auto& test : test_cases) {
        if (test.first.length() == test.second.length()) {
            double identity = calculateIdentity(test.first, test.second);
            EXPECT_GE(identity, 0.0);
            EXPECT_LE(identity, 100.0);
        }
    }
}

