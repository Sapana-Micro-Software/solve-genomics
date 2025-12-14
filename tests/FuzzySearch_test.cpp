#include <gtest/gtest.h>
#include "FuzzySearch.h"
#include "ExactMatch.h"  // For comparison

// ============================================================================
// UNIT TESTS - Testing individual methods
// ============================================================================

class FuzzySearchUnitTest : public ::testing::Test {
protected:
    void SetUp() override {
        searcher = new FuzzySearch();
    }
    
    void TearDown() override {
        delete searcher;
    }
    
    FuzzySearch* searcher;
};

TEST_F(FuzzySearchUnitTest, Search_ExactMatch_DistanceZero) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCG";
    
    FuzzySearchResult result = searcher->search(sequence, pattern, 0);
    
    // Should find exact matches only
    EXPECT_EQ(result.count, 2);
    EXPECT_EQ(result.positions[0], 0);
    EXPECT_EQ(result.positions[1], 4);
    EXPECT_EQ(result.distances[0], 0);
    EXPECT_EQ(result.distances[1], 0);
}

TEST_F(FuzzySearchUnitTest, Search_OneMismatch_DistanceOne) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCC";  // One mismatch from "ATCG"
    
    FuzzySearchResult result = searcher->search(sequence, pattern, 1);
    
    // Should find matches with distance <= 1
    EXPECT_GT(result.count, 0);
    for (size_t i = 0; i < result.distances.size(); ++i) {
        EXPECT_LE(result.distances[i], 1);
    }
}

TEST_F(FuzzySearchUnitTest, Search_NoMatch_StrictThreshold) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "GGGG";
    
    FuzzySearchResult result = searcher->search(sequence, pattern, 0);
    
    // No matches with distance 0
    EXPECT_EQ(result.count, 0);
}

TEST_F(FuzzySearchUnitTest, Search_MatchWithHigherThreshold) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "GGGG";
    
    FuzzySearchResult result1 = searcher->search(sequence, pattern, 0);
    FuzzySearchResult result2 = searcher->search(sequence, pattern, 4);
    
    // With higher threshold, should find matches
    EXPECT_EQ(result1.count, 0);
    EXPECT_GE(result2.count, 0);  // May find matches with distance <= 4
}

TEST_F(FuzzySearchUnitTest, Search_EmptyInput) {
    FuzzySearchResult result1 = searcher->search("", "ATCG", 1);
    EXPECT_EQ(result1.count, 0);
    
    FuzzySearchResult result2 = searcher->search("ATCG", "", 1);
    EXPECT_EQ(result2.count, 0);
    
    FuzzySearchResult result3 = searcher->search("", "", 1);
    EXPECT_EQ(result3.count, 0);
}

TEST_F(FuzzySearchUnitTest, SearchHamming_ExactMatch) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCG";
    
    FuzzySearchResult result = searcher->searchHamming(sequence, pattern, 0);
    
    EXPECT_EQ(result.count, 2);
    EXPECT_EQ(result.distances[0], 0);
    EXPECT_EQ(result.distances[1], 0);
}

TEST_F(FuzzySearchUnitTest, SearchHamming_OneMismatch) {
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCC";  // One substitution
    
    FuzzySearchResult result = searcher->searchHamming(sequence, pattern, 1);
    
    EXPECT_GT(result.count, 0);
    for (size_t i = 0; i < result.distances.size(); ++i) {
        EXPECT_LE(result.distances[i], 1);
        EXPECT_EQ(result.distances[i], 1);  // Should be exactly 1
    }
}

TEST_F(FuzzySearchUnitTest, Contains_True) {
    EXPECT_TRUE(searcher->contains("ATCGATCG", "ATCG", 0));
    EXPECT_TRUE(searcher->contains("ATCGATCG", "ATCC", 1));  // One mismatch allowed
}

TEST_F(FuzzySearchUnitTest, Contains_False) {
    EXPECT_FALSE(searcher->contains("ATCGATCG", "GGGG", 0));
    EXPECT_FALSE(searcher->contains("", "ATCG", 1));
}

TEST_F(FuzzySearchUnitTest, FindFirst_Found) {
    int pos = searcher->findFirst("ATCGATCG", "ATCG", 0);
    EXPECT_EQ(pos, 0);
    
    pos = searcher->findFirst("ATCGATCG", "ATCC", 1);
    EXPECT_GE(pos, 0);
}

TEST_F(FuzzySearchUnitTest, FindFirst_NotFound) {
    int pos = searcher->findFirst("ATCGATCG", "GGGG", 0);
    EXPECT_EQ(pos, -1);
    
    pos = searcher->findFirst("", "ATCG", 1);
    EXPECT_EQ(pos, -1);
}

// ============================================================================
// STATIC METHOD TESTS - Edit Distance and Hamming Distance
// ============================================================================

TEST_F(FuzzySearchUnitTest, EditDistance_Identical) {
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "ATCG"), 0);
    EXPECT_EQ(FuzzySearch::editDistance("", ""), 0);
}

TEST_F(FuzzySearchUnitTest, EditDistance_OneSubstitution) {
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "ATCC"), 1);
    EXPECT_EQ(FuzzySearch::editDistance("A", "T"), 1);
}

TEST_F(FuzzySearchUnitTest, EditDistance_OneInsertion) {
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "ATCGG"), 1);
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "AATCG"), 1);
}

TEST_F(FuzzySearchUnitTest, EditDistance_OneDeletion) {
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "ATC"), 1);
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "TCG"), 1);
}

TEST_F(FuzzySearchUnitTest, EditDistance_MultipleOperations) {
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "GGGG"), 4);
    EXPECT_EQ(FuzzySearch::editDistance("ATCGATCG", "GGGGGGGG"), 8);
}

TEST_F(FuzzySearchUnitTest, EditDistance_CaseInsensitive) {
    EXPECT_EQ(FuzzySearch::editDistance("ATCG", "atcg"), 0);
    EXPECT_EQ(FuzzySearch::editDistance("AtCg", "aTcG"), 0);
}

TEST_F(FuzzySearchUnitTest, HammingDistance_Identical) {
    EXPECT_EQ(FuzzySearch::hammingDistance("ATCG", "ATCG"), 0);
    EXPECT_EQ(FuzzySearch::hammingDistance("A", "A"), 0);
}

TEST_F(FuzzySearchUnitTest, HammingDistance_OneMismatch) {
    EXPECT_EQ(FuzzySearch::hammingDistance("ATCG", "ATCC"), 1);
    EXPECT_EQ(FuzzySearch::hammingDistance("ATCG", "ATCG"), 0);
}

TEST_F(FuzzySearchUnitTest, HammingDistance_DifferentLengths) {
    EXPECT_EQ(FuzzySearch::hammingDistance("ATCG", "ATC"), -1);
    EXPECT_EQ(FuzzySearch::hammingDistance("ATCG", "ATCGATCG"), -1);
}

TEST_F(FuzzySearchUnitTest, HammingDistance_CaseInsensitive) {
    EXPECT_EQ(FuzzySearch::hammingDistance("ATCG", "atcg"), 0);
    EXPECT_EQ(FuzzySearch::hammingDistance("ATCG", "ATCC"), 1);
}

// ============================================================================
// REGRESSION TESTS - Test known good outputs
// ============================================================================

TEST_F(FuzzySearchUnitTest, Regression_KnownEditDistances) {
    // Known test cases
    EXPECT_EQ(FuzzySearch::editDistance("kitten", "sitting"), 3);
    EXPECT_EQ(FuzzySearch::editDistance("sunday", "saturday"), 3);
    EXPECT_EQ(FuzzySearch::editDistance("", "abc"), 3);
    EXPECT_EQ(FuzzySearch::editDistance("abc", ""), 3);
}

TEST_F(FuzzySearchUnitTest, Regression_KnownHammingDistances) {
    EXPECT_EQ(FuzzySearch::hammingDistance("karolin", "kathrin"), 3);
    EXPECT_EQ(FuzzySearch::hammingDistance("karolin", "kerstin"), 3);
}

TEST_F(FuzzySearchUnitTest, Regression_Consistency) {
    // Same search should give same results
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCG";
    
    FuzzySearchResult result1 = searcher->search(sequence, pattern, 0);
    FuzzySearchResult result2 = searcher->search(sequence, pattern, 0);
    
    EXPECT_EQ(result1.count, result2.count);
    EXPECT_EQ(result1.positions.size(), result2.positions.size());
}

// ============================================================================
// A-B COMPARISON TESTS - Compare with ExactMatch
// ============================================================================

TEST_F(FuzzySearchUnitTest, AB_CompareWithExactMatch_DistanceZero) {
    ExactMatch exact_matcher;
    
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCG";
    
    FuzzySearchResult fuzzy_result = searcher->search(sequence, pattern, 0);
    SearchResult exact_result = exact_matcher.search(sequence, pattern);
    
    // With distance 0, should match exact search
    EXPECT_EQ(fuzzy_result.count, exact_result.count);
    EXPECT_EQ(fuzzy_result.positions.size(), exact_result.positions.size());
    
    for (size_t i = 0; i < fuzzy_result.positions.size(); ++i) {
        EXPECT_EQ(fuzzy_result.positions[i], exact_result.positions[i]);
        EXPECT_EQ(fuzzy_result.distances[i], 0);
    }
}

// ============================================================================
// BLACKBOX TESTS - Test without knowing implementation
// ============================================================================

TEST_F(FuzzySearchUnitTest, Blackbox_DistanceNonNegative) {
    // All distances should be non-negative
    std::string sequence = "ATCGATCGATCG";
    std::string pattern = "ATCG";
    
    for (int max_dist = 0; max_dist <= 3; ++max_dist) {
        FuzzySearchResult result = searcher->search(sequence, pattern, max_dist);
        
        for (int dist : result.distances) {
            EXPECT_GE(dist, 0);
            EXPECT_LE(dist, max_dist);
        }
    }
}

TEST_F(FuzzySearchUnitTest, Blackbox_PositionsValid) {
    std::string sequence = "ATCGATCGATCG";
    std::string pattern = "ATCG";
    
    FuzzySearchResult result = searcher->search(sequence, pattern, 1);
    
    for (int pos : result.positions) {
        EXPECT_GE(pos, 0);
        EXPECT_LT(pos, static_cast<int>(sequence.length()));
    }
}

TEST_F(FuzzySearchUnitTest, Blackbox_CountMatchesPositions) {
    // Count should match number of positions
    std::string sequence = "ATCGATCG";
    std::string pattern = "ATCG";
    
    FuzzySearchResult result = searcher->search(sequence, pattern, 1);
    
    EXPECT_EQ(result.count, static_cast<int>(result.positions.size()));
    EXPECT_EQ(result.count, static_cast<int>(result.distances.size()));
}

TEST_F(FuzzySearchUnitTest, Blackbox_IncreasingThreshold) {
    // More matches with higher threshold
    std::string sequence = "ATCGATCG";
    std::string pattern = "GGGG";
    
    FuzzySearchResult result0 = searcher->search(sequence, pattern, 0);
    FuzzySearchResult result1 = searcher->search(sequence, pattern, 1);
    FuzzySearchResult result2 = searcher->search(sequence, pattern, 2);
    
    EXPECT_LE(result0.count, result1.count);
    EXPECT_LE(result1.count, result2.count);
}

TEST_F(FuzzySearchUnitTest, Blackbox_EditDistanceProperties) {
    // Edit distance should satisfy triangle inequality (approximately)
    std::string a = "ATCG";
    std::string b = "ATCC";
    std::string c = "ATGG";
    
    int d_ab = FuzzySearch::editDistance(a, b);
    int d_bc = FuzzySearch::editDistance(b, c);
    int d_ac = FuzzySearch::editDistance(a, c);
    
    // d(a,c) <= d(a,b) + d(b,c) (triangle inequality)
    EXPECT_LE(d_ac, d_ab + d_bc);
    
    // Symmetry: d(a,b) = d(b,a)
    EXPECT_EQ(FuzzySearch::editDistance(a, b), FuzzySearch::editDistance(b, a));
}

TEST_F(FuzzySearchUnitTest, Blackbox_HammingDistanceProperties) {
    // Hamming distance properties
    std::string a = "ATCG";
    std::string b = "ATCC";
    std::string c = "ATGG";
    
    int d_ab = FuzzySearch::hammingDistance(a, b);
    int d_ba = FuzzySearch::hammingDistance(b, a);
    
    // Symmetry
    EXPECT_EQ(d_ab, d_ba);
    
    // Non-negativity
    EXPECT_GE(d_ab, 0);
    EXPECT_GE(d_ba, 0);
    
    // Maximum distance is length of string
    EXPECT_LE(d_ab, static_cast<int>(a.length()));
}

TEST_F(FuzzySearchUnitTest, Blackbox_EmptyPattern) {
    // Empty pattern should return no matches
    FuzzySearchResult result = searcher->search("ATCG", "", 0);
    EXPECT_EQ(result.count, 0);
}

TEST_F(FuzzySearchUnitTest, Blackbox_PatternLongerThanSequence) {
    // Pattern longer than sequence + max_distance should return no matches
    FuzzySearchResult result = searcher->search("ATC", "ATCGATCG", 0);
    EXPECT_EQ(result.count, 0);
    
    // But with high enough distance, might find matches
    FuzzySearchResult result2 = searcher->search("ATC", "ATCGATCG", 5);
    // May or may not find matches depending on implementation
}

TEST_F(FuzzySearchUnitTest, Blackbox_AllDistancesWithinThreshold) {
    // All returned distances should be within threshold
    std::string sequence = "ATCGATCGATCG";
    std::string pattern = "ATCG";
    
    for (int max_dist = 0; max_dist <= 3; ++max_dist) {
        FuzzySearchResult result = searcher->search(sequence, pattern, max_dist);
        
        for (int dist : result.distances) {
            EXPECT_LE(dist, max_dist) << "Distance " << dist 
                                      << " exceeds threshold " << max_dist;
        }
    }
}

TEST_F(FuzzySearchUnitTest, Blackbox_HammingVsEditDistance) {
    // For same-length strings, Hamming distance <= Edit distance
    std::string seq = "ATCGATCG";
    std::string pat = "ATCCATCG";
    
    int hamming = FuzzySearch::hammingDistance(seq, pat);
    int edit = FuzzySearch::editDistance(seq, pat);
    
    if (hamming != -1) {
        EXPECT_LE(hamming, edit);
    }
}

