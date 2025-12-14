#include "association_list.h"
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <vector>

void AssociationList::addAssociation(const std::string& symbol, const std::string& sequence) {
    associations_[symbol] = sequence;
}

std::string AssociationList::getSequence(const std::string& symbol) const {
    auto it = associations_.find(symbol);
    if (it != associations_.end()) {
        return it->second;
    }
    return "";  // Symbol not found
}

bool AssociationList::hasSymbol(const std::string& symbol) const {
    return associations_.find(symbol) != associations_.end();
}

std::string AssociationList::expand(const std::string& compressed_sequence) const {
    std::string expanded = compressed_sequence;
    
    // Replace symbols with their sequences (process longest symbols first)
    std::vector<std::pair<std::string, std::string>> sorted_associations;
    for (const auto& assoc : associations_) {
        sorted_associations.push_back(assoc);
    }
    
    std::sort(sorted_associations.begin(), sorted_associations.end(),
        [](const std::pair<std::string, std::string>& a,
           const std::pair<std::string, std::string>& b) {
            return a.first.length() > b.first.length();
        });
    
    for (const auto& assoc : sorted_associations) {
        const std::string& symbol = assoc.first;
        const std::string& sequence = assoc.second;
        
        size_t pos = 0;
        while ((pos = expanded.find(symbol, pos)) != std::string::npos) {
            expanded.replace(pos, symbol.length(), sequence);
            pos += sequence.length();
        }
    }
    
    return expanded;
}

void AssociationList::buildFromGrammar(const std::vector<std::string>& grammar) {
    associations_.clear();
    
    for (const std::string& rule : grammar) {
        size_t arrow_pos = rule.find(" -> ");
        if (arrow_pos != std::string::npos) {
            std::string symbol = rule.substr(0, arrow_pos);
            std::string sequence = rule.substr(arrow_pos + 4);
            associations_[symbol] = sequence;
        }
    }
}

void AssociationList::buildFromDictionary(const std::map<std::string, std::string>& dictionary) {
    associations_ = dictionary;
}

std::vector<std::string> AssociationList::getSymbols() const {
    std::vector<std::string> symbols;
    for (const auto& assoc : associations_) {
        symbols.push_back(assoc.first);
    }
    return symbols;
}

void AssociationList::clear() {
    associations_.clear();
}

size_t AssociationList::size() const {
    return associations_.size();
}

// CompressedPatternMatcher implementation

CompressedPatternMatcher::CompressedPatternMatcher(const AssociationList& association_list)
    : association_list_(association_list) {
}

std::vector<size_t> CompressedPatternMatcher::search(
    const std::string& compressed_sequence,
    const std::string& pattern) {
    
    // Expand compressed sequence
    std::string expanded = expandSequence(compressed_sequence);
    
    // Search in expanded sequence
    std::vector<size_t> positions;
    size_t pos = 0;
    while ((pos = expanded.find(pattern, pos)) != std::string::npos) {
        positions.push_back(pos);
        pos += pattern.length();
    }
    
    return positions;
}

std::vector<std::pair<size_t, int>> CompressedPatternMatcher::searchFuzzy(
    const std::string& compressed_sequence,
    const std::string& pattern,
    int max_distance) {
    
    // For fuzzy search, we need to expand and use edit distance
    std::string expanded = expandSequence(compressed_sequence);
    
    std::vector<std::pair<size_t, int>> results;
    
    // Slide pattern over expanded sequence
    for (size_t i = 0; i <= expanded.length() - pattern.length() + max_distance; ++i) {
        size_t end_pos = std::min(i + pattern.length() + max_distance, expanded.length());
        std::string substring = expanded.substr(i, end_pos - i);
        
        // Calculate edit distance
        int distance = calculateEditDistance(pattern, substring, max_distance);
        
        if (distance <= max_distance) {
            results.push_back({i, distance});
        }
    }
    
    return results;
}

std::string CompressedPatternMatcher::expandSequence(const std::string& compressed) const {
    return association_list_.expand(compressed);
}

size_t CompressedPatternMatcher::convertPosition(
    const std::string& compressed,
    size_t compressed_pos) const {
    
    // Convert position in compressed sequence to original sequence
    std::string expanded = expandSequence(compressed);
    // For simplicity, return the position (would need more complex mapping)
    return compressed_pos;
}

// Helper function for edit distance (simplified version)
int CompressedPatternMatcher::calculateEditDistance(
    const std::string& str1,
    const std::string& str2,
    int max_distance) const {
    
    int m = str1.length();
    int n = str2.length();
    
    // Early termination
    int len_diff = (m > n) ? (m - n) : (n - m);
    if (len_diff > max_distance) {
        return max_distance + 1;
    }
    
    // Simplified DP (could use FuzzySearch::editDistanceBounded)
    std::vector<int> prev(n + 1);
    std::vector<int> curr(n + 1);
    
    for (int j = 0; j <= n; ++j) {
        prev[j] = j;
    }
    
    for (int i = 1; i <= m; ++i) {
        curr[0] = i;
        for (int j = 1; j <= n; ++j) {
            if (std::toupper(str1[i - 1]) == std::toupper(str2[j - 1])) {
                curr[j] = prev[j - 1];
            } else {
                curr[j] = 1 + std::min({prev[j], curr[j - 1], prev[j - 1]});
            }
        }
        
        int min_in_row = *std::min_element(curr.begin(), curr.end());
        if (min_in_row > max_distance) {
            return max_distance + 1;
        }
        
        prev.swap(curr);
    }
    
    return prev[n];
}

