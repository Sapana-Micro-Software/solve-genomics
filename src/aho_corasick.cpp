#include "aho_corasick.h"
#include <algorithm>
#include <cctype>

AhoCorasick::AhoCorasick() : node_count_(0) {
    trie_.push_back(TrieNode());  // Root node
    node_count_ = 1;
}

void AhoCorasick::buildAutomaton(const std::vector<std::string>& patterns) {
    // Clear previous automaton
    trie_.clear();
    trie_.push_back(TrieNode());
    node_count_ = 1;
    patterns_ = patterns;
    
    // Add all patterns to trie
    for (size_t i = 0; i < patterns.size(); ++i) {
        addPattern(patterns[i], static_cast<int>(i));
    }
    
    // Build failure links
    buildFailureLinks();
}

void AhoCorasick::addPattern(const std::string& pattern, int pattern_index) {
    int current = 0;  // Root
    
    for (char c : pattern) {
        char upper_c = std::toupper(c);
        
        if (trie_[current].children.find(upper_c) == trie_[current].children.end()) {
            // Create new node
            trie_.push_back(TrieNode());
            trie_[current].children[upper_c] = node_count_;
            current = node_count_;
            node_count_++;
        } else {
            current = trie_[current].children[upper_c];
        }
    }
    
    // Mark end of pattern
    trie_[current].is_end = true;
    trie_[current].output.push_back(pattern_index);
}

void AhoCorasick::buildFailureLinks() {
    std::queue<int> q;
    
    // Initialize root's children failure links to root
    for (const auto& child_pair : trie_[0].children) {
        trie_[child_pair.second].fail = 0;
        q.push(child_pair.second);
    }
    
    // BFS to build failure links
    while (!q.empty()) {
        int current = q.front();
        q.pop();
        
        // For each child of current node
        for (const auto& child_pair : trie_[current].children) {
            char c = child_pair.first;
            int child = child_pair.second;
            // Find failure link for child
            int fail = trie_[current].fail;
            
            while (fail != -1 && trie_[fail].children.find(c) == trie_[fail].children.end()) {
                fail = trie_[fail].fail;
            }
            
            if (fail == -1) {
                trie_[child].fail = 0;  // Root
            } else {
                trie_[child].fail = trie_[fail].children[c];
            }
            
            // Merge output from failure link
            if (trie_[child].fail != -1) {
                trie_[child].output.insert(trie_[child].output.end(),
                                          trie_[trie_[child].fail].output.begin(),
                                          trie_[trie_[child].fail].output.end());
            }
            
            q.push(child);
        }
    }
}

AhoCorasick::MultiPatternResult AhoCorasick::search(const std::string& text) {
    MultiPatternResult result;
    
    if (trie_.empty() || text.empty()) {
        return result;
    }
    
    int current = 0;  // Root
    
    for (size_t i = 0; i < text.length(); ++i) {
        char c = std::toupper(text[i]);
        
        // Follow failure links until we find a valid transition
        while (current != -1 && trie_[current].children.find(c) == trie_[current].children.end()) {
            current = trie_[current].fail;
        }
        
        if (current == -1) {
            current = 0;  // Root
        } else {
            current = trie_[current].children[c];
        }
        
        // Collect all matches at current node
        for (int pattern_idx : trie_[current].output) {
            size_t pattern_len = patterns_[pattern_idx].length();
            size_t match_pos = i - pattern_len + 1;
            
            result.matches[patterns_[pattern_idx]].push_back(match_pos);
            result.total_matches++;
        }
    }
    
    return result;
}

SearchResult AhoCorasick::searchSingle(const std::string& text, const std::string& pattern) {
    SearchResult result;
    
    // Build automaton with single pattern
    buildAutomaton({pattern});
    
    // Search
    MultiPatternResult multi_result = search(text);
    
    if (multi_result.matches.find(pattern) != multi_result.matches.end()) {
        result.positions.clear();
        for (size_t pos : multi_result.matches[pattern]) {
            result.positions.push_back(static_cast<int>(pos));
        }
        result.count = static_cast<int>(result.positions.size());
    }
    
    return result;
}

int AhoCorasick::getFailureLink(int node, char c) {
    while (node != -1 && trie_[node].children.find(c) == trie_[node].children.end()) {
        node = trie_[node].fail;
    }
    
    if (node == -1) {
        return 0;  // Root
    }
    
    return trie_[node].children[c];
}

void AhoCorasick::collectOutputs(int node, std::vector<int>& outputs) {
    outputs.insert(outputs.end(), trie_[node].output.begin(), trie_[node].output.end());
}

