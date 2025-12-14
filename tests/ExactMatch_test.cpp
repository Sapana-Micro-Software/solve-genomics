#include <gtest/gtest.h>
#include "ExactMatch.h"
#include "utils.h"

// ============================================================================
// UNIT TESTS - Testing individual methods
// ============================================================================

class ExactMatchUnitTest : public ::testing::Test {
protected:
    void SetUp() override {
        matcher = new ExactMatch();
    }
    
    void TearDown() override {
        delete matcher;
    }
    
    ExactMatch* matcher;
};

TEST_F(ExactMatchUnitTest, Search_SingleMatch) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCG";
    
    SearchResult result = matcher->search(sequence, pattern);
    
    EXPECT_EQ(result.count, 2);
    EXPECT_EQ(result.positions.size(), 2);
    EXPECT_EQ(result.positions[0], 0);
    EXPECT_EQ(result.positions[1], 4);
}

TEST_F(ExactMatchUnitTest, Search_NoMatch) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "GGGG";
    
    SearchResult result = matcher->search(sequence, pattern);
    
    EXPECT_EQ(result.count, 0);
    EXPECT_TRUE(result.positions.empty());
}

TEST_F(ExactMatchUnitTest, Search_EmptyInput) {
    SearchResult result1 = matcher->search("", "ATCG");
    EXPECT_EQ(result1.count, 0);
    
    SearchResult result2 = matcher->search("ATCG", "");
    EXPECT_EQ(result2.count, 0);
    
    SearchResult result3 = matcher->search("", "");
    EXPECT_EQ(result3.count, 0);
}

TEST_F(ExactMatchUnitTest, Search_PatternLongerThanSequence) {
    std::string sequence = "ATC";
    std::string pattern = "ATCGATCG";
    
    SearchResult result = matcher->search(sequence, pattern);
    EXPECT_EQ(result.count, 0);
}

TEST_F(ExactMatchUnitTest, Search_CaseInsensitive) {
    std::string sequence = "atcgatcg";
    std::string pattern = "ATCG";
    
    SearchResult result = matcher->search(sequence, pattern);
    EXPECT_GT(result.count, 0);
}

TEST_F(ExactMatchUnitTest, Contains_True) {
    EXPECT_TRUE(matcher->contains("ATCGATCG", "ATCG"));
    EXPECT_TRUE(matcher->contains("ATCGATCG", "TCGA"));
}

TEST_F(ExactMatchUnitTest, Contains_False) {
    EXPECT_FALSE(matcher->contains("ATCGATCG", "GGGG"));
    EXPECT_FALSE(matcher->contains("", "ATCG"));
}

TEST_F(ExactMatchUnitTest, FindFirst_Found) {
    int pos = matcher->findFirst("ATCGATCG", "ATCG");
    EXPECT_EQ(pos, 0);
    
    pos = matcher->findFirst("ATCGATCG", "TCGA");
    EXPECT_EQ(pos, 1);
}

TEST_F(ExactMatchUnitTest, FindFirst_NotFound) {
    int pos = matcher->findFirst("ATCGATCG", "GGGG");
    EXPECT_EQ(pos, -1);
    
    pos = matcher->findFirst("", "ATCG");
    EXPECT_EQ(pos, -1);
}

// ============================================================================
// REGRESSION TESTS - Test known good outputs
// ============================================================================

TEST_F(ExactMatchUnitTest, Regression_KnownSequences) {
    // Known test case: pattern appears multiple times
    std::string sequence = "ATCGATCGATCGATCG";
    std::string pattern = "ATCG";
    
    SearchResult result = matcher->search(sequence, pattern);
    
    // Should find 4 occurrences at positions 0, 4, 8, 12
    EXPECT_EQ(result.count, 4);
    EXPECT_EQ(result.positions[0], 0);
    EXPECT_EQ(result.positions[1], 4);
    EXPECT_EQ(result.positions[2], 8);
    EXPECT_EQ(result.positions[3], 12);
}

TEST_F(ExactMatchUnitTest, Regression_OverlappingPatterns) {
    // Pattern "AAA" in "AAAAA" should find 3 matches (positions 0, 1, 2)
    std::string sequence = "AAAAA";
    std::string pattern = "AAA";
    
    SearchResult result = matcher->search(sequence, pattern);
    EXPECT_EQ(result.count, 3);
    EXPECT_EQ(result.positions[0], 0);
    EXPECT_EQ(result.positions[1], 1);
    EXPECT_EQ(result.positions[2], 2);
}

// ============================================================================
// BLACKBOX TESTS - Test without knowing implementation
// ============================================================================

TEST_F(ExactMatchUnitTest, Blackbox_Consistency) {
    // Search and contains should be consistent
    std::string sequence = "ATCGATCGATCG";
    std::string pattern = "ATCG";
    
    bool contains_result = matcher->contains(sequence, pattern);
    SearchResult search_result = matcher->search(sequence, pattern);
    int find_result = matcher->findFirst(sequence, pattern);
    
    EXPECT_EQ(contains_result, search_result.count > 0);
    if (search_result.count > 0) {
        EXPECT_EQ(find_result, search_result.positions[0]);
    } else {
        EXPECT_EQ(find_result, -1);
    }
}

TEST_F(ExactMatchUnitTest, Blackbox_AllPositionsValid) {
    // All returned positions should be valid indices
    std::string sequence = "ATCGATCGATCGATCGATCG";
    std::string pattern = "ATCG";
    
    SearchResult result = matcher->search(sequence, pattern);
    
    for (int pos : result.positions) {
        EXPECT_GE(pos, 0);
        EXPECT_LT(pos, static_cast<int>(sequence.length()));
        
        // Verify the pattern actually matches at this position
        std::string substring = sequence.substr(pos, pattern.length());
        EXPECT_EQ(substring.length(), pattern.length());
        // Case-insensitive comparison
        for (size_t i = 0; i < pattern.length(); ++i) {
            EXPECT_EQ(std::toupper(substring[i]), std::toupper(pattern[i]));
        }
    }
}

TEST_F(ExactMatchUnitTest, Blackbox_CountMatchesPositions) {
    // Count should match the number of positions
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATCGATCG", "ATCG"},
        {"AAAAA", "AAA"},
        {"ATCG", "GGGG"},
        {"A", "A"},
        {"ATCGATCG", "ATCGATCG"}
    };
    
    for (const auto& test : test_cases) {
        SearchResult result = matcher->search(test.first, test.second);
        EXPECT_EQ(result.count, static_cast<int>(result.positions.size()));
    }
}

TEST_F(ExactMatchUnitTest, Blackbox_EmptyPattern) {
    // Empty pattern should return no matches
    SearchResult result = matcher->search("ATCG", "");
    EXPECT_EQ(result.count, 0);
    EXPECT_TRUE(result.positions.empty());
}

TEST_F(ExactMatchUnitTest, Blackbox_ExactMatchProperty) {
    // If pattern found at position i, then sequence[i:i+len(pattern)] should equal pattern
    std::string sequence = "ATCGATCGATCG";
    std::string pattern = "TCGA";
    
    SearchResult result = matcher->search(sequence, pattern);
    
    for (int pos : result.positions) {
        std::string substring = sequence.substr(pos, pattern.length());
        // Case-insensitive comparison
        bool matches = true;
        for (size_t i = 0; i < pattern.length(); ++i) {
            if (std::toupper(substring[i]) != std::toupper(pattern[i])) {
                matches = false;
                break;
            }
        }
        EXPECT_TRUE(matches);
    }
}

