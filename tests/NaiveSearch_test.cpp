#include <gtest/gtest.h>
#include "NaiveSearch.h"
#include "ExactMatch.h"  // For comparison

// ============================================================================
// UNIT TESTS - Testing individual methods
// ============================================================================

class NaiveSearchUnitTest : public ::testing::Test {
protected:
    void SetUp() override {
        searcher = new NaiveSearch();
    }
    
    void TearDown() override {
        delete searcher;
    }
    
    NaiveSearch* searcher;
};

TEST_F(NaiveSearchUnitTest, Search_SingleMatch) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCG";
    
    SearchResult result = searcher->search(sequence, pattern);
    
    EXPECT_EQ(result.count, 2);
    EXPECT_EQ(result.positions.size(), 2);
    EXPECT_EQ(result.positions[0], 0);
    EXPECT_EQ(result.positions[1], 4);
}

TEST_F(NaiveSearchUnitTest, Search_NoMatch) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "GGGG";
    
    SearchResult result = searcher->search(sequence, pattern);
    
    EXPECT_EQ(result.count, 0);
    EXPECT_TRUE(result.positions.empty());
}

TEST_F(NaiveSearchUnitTest, Search_EmptyInput) {
    SearchResult result1 = searcher->search("", "ATCG");
    EXPECT_EQ(result1.count, 0);
    
    SearchResult result2 = searcher->search("ATCG", "");
    EXPECT_EQ(result2.count, 0);
    
    SearchResult result3 = searcher->search("", "");
    EXPECT_EQ(result3.count, 0);
}

TEST_F(NaiveSearchUnitTest, Search_PatternLongerThanSequence) {
    std::string sequence = "ATC";
    std::string pattern = "ATCGATCG";
    
    SearchResult result = searcher->search(sequence, pattern);
    EXPECT_EQ(result.count, 0);
}

TEST_F(NaiveSearchUnitTest, Search_CaseInsensitive) {
    std::string sequence = "atcgatcg";
    std::string pattern = "ATCG";
    
    SearchResult result = searcher->search(sequence, pattern);
    EXPECT_GT(result.count, 0);
}

TEST_F(NaiveSearchUnitTest, FindFirst_Found) {
    int pos = searcher->findFirst("ATCGATCG", "ATCG");
    EXPECT_EQ(pos, 0);
    
    pos = searcher->findFirst("ATCGATCG", "TCGA");
    EXPECT_EQ(pos, 1);
}

TEST_F(NaiveSearchUnitTest, FindFirst_NotFound) {
    int pos = searcher->findFirst("ATCGATCG", "GGGG");
    EXPECT_EQ(pos, -1);
    
    pos = searcher->findFirst("", "ATCG");
    EXPECT_EQ(pos, -1);
}

// ============================================================================
// REGRESSION TESTS - Test known good outputs
// ============================================================================

TEST_F(NaiveSearchUnitTest, Regression_KnownSequences) {
    std::string sequence = "ATCGATCGATCGATCG";
    std::string pattern = "ATCG";
    
    SearchResult result = searcher->search(sequence, pattern);
    
    EXPECT_EQ(result.count, 4);
    EXPECT_EQ(result.positions[0], 0);
    EXPECT_EQ(result.positions[1], 4);
    EXPECT_EQ(result.positions[2], 8);
    EXPECT_EQ(result.positions[3], 12);
}

TEST_F(NaiveSearchUnitTest, Regression_OverlappingPatterns) {
    std::string sequence = "AAAAA";
    std::string pattern = "AAA";
    
    SearchResult result = searcher->search(sequence, pattern);
    EXPECT_EQ(result.count, 3);
    EXPECT_EQ(result.positions[0], 0);
    EXPECT_EQ(result.positions[1], 1);
    EXPECT_EQ(result.positions[2], 2);
}

// ============================================================================
// A-B TESTS - Compare with ExactMatch (should produce same results)
// ============================================================================

TEST_F(NaiveSearchUnitTest, AB_CompareWithExactMatch) {
    ExactMatch exact_matcher;
    
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATCGATCG", "ATCG"},
        {"AAAAA", "AAA"},
        {"ATCG", "GGGG"},
        {"ATCGATCGATCG", "TCGA"},
        {"A", "A"},
        {"ATCGATCG", "ATCGATCG"}
    };
    
    for (const auto& test : test_cases) {
        SearchResult naive_result = searcher->search(test.first, test.second);
        SearchResult exact_result = exact_matcher.search(test.first, test.second);
        
        EXPECT_EQ(naive_result.count, exact_result.count);
        EXPECT_EQ(naive_result.positions.size(), exact_result.positions.size());
        
        for (size_t i = 0; i < naive_result.positions.size(); ++i) {
            EXPECT_EQ(naive_result.positions[i], exact_result.positions[i]);
        }
    }
}

// ============================================================================
// BLACKBOX TESTS - Test without knowing implementation
// ============================================================================

TEST_F(NaiveSearchUnitTest, Blackbox_AllPositionsValid) {
    std::string sequence = "ATCGATCGATCGATCGATCG";
    std::string pattern = "ATCG";
    
    SearchResult result = searcher->search(sequence, pattern);
    
    for (int pos : result.positions) {
        EXPECT_GE(pos, 0);
        EXPECT_LT(pos, static_cast<int>(sequence.length()));
        
        std::string substring = sequence.substr(pos, pattern.length());
        EXPECT_EQ(substring.length(), pattern.length());
        for (size_t i = 0; i < pattern.length(); ++i) {
            EXPECT_EQ(std::toupper(substring[i]), std::toupper(pattern[i]));
        }
    }
}

TEST_F(NaiveSearchUnitTest, Blackbox_CountMatchesPositions) {
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATCGATCG", "ATCG"},
        {"AAAAA", "AAA"},
        {"ATCG", "GGGG"},
        {"A", "A"}
    };
    
    for (const auto& test : test_cases) {
        SearchResult result = searcher->search(test.first, test.second);
        EXPECT_EQ(result.count, static_cast<int>(result.positions.size()));
    }
}

TEST_F(NaiveSearchUnitTest, Blackbox_ExactMatchProperty) {
    std::string sequence = "ATCGATCGATCG";
    std::string pattern = "TCGA";
    
    SearchResult result = searcher->search(sequence, pattern);
    
    for (int pos : result.positions) {
        std::string substring = sequence.substr(pos, pattern.length());
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

TEST_F(NaiveSearchUnitTest, Blackbox_NoFalsePositives) {
    // Verify that returned positions actually contain the pattern
    std::string sequence = "ATCGATCGATCG";
    std::string pattern = "ATCG";
    
    SearchResult result = searcher->search(sequence, pattern);
    
    for (int pos : result.positions) {
        std::string substring = sequence.substr(pos, pattern.length());
        EXPECT_EQ(substring.length(), pattern.length());
        
        // Check character by character (case-insensitive)
        for (size_t i = 0; i < pattern.length(); ++i) {
            EXPECT_EQ(std::toupper(substring[i]), std::toupper(pattern[i]))
                << "Position " << pos << " does not match pattern at index " << i;
        }
    }
}

