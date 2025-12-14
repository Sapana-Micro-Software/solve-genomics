#ifndef AHO_CORASICK_H
#define AHO_CORASICK_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>
#include <queue>

/**
 * Aho-Corasick Algorithm for Multi-Pattern Matching
 * Builds a finite automaton to search for multiple patterns simultaneously
 * Time: O(n + m + z) where n is text length, m is total pattern length, z is matches
 * Space: O(m)
 */
class AhoCorasick {
public:
    struct MultiPatternResult {
        std::map<std::string, std::vector<size_t>> matches;  // Pattern -> positions
        int total_matches;
        
        MultiPatternResult() : total_matches(0) {}
    };
    
    AhoCorasick();
    
    /**
     * Build automaton from multiple patterns
     * @param patterns Vector of patterns to search for
     */
    void buildAutomaton(const std::vector<std::string>& patterns);
    
    /**
     * Search for all patterns in text
     * @param text Text to search in
     * @return MultiPatternResult with all matches
     */
    MultiPatternResult search(const std::string& text);
    
    /**
     * Search for single pattern (for compatibility)
     * @param text Text to search in
     * @param pattern Pattern to find
     * @return SearchResult with matches
     */
    SearchResult searchSingle(const std::string& text, const std::string& pattern);
    
private:
    struct TrieNode {
        std::map<char, int> children;  // Character -> node index
        std::vector<int> output;       // Pattern indices ending at this node
        int fail;                      // Failure link
        bool is_end;                   // End of pattern flag
        
        TrieNode() : fail(-1), is_end(false) {}
    };
    
    std::vector<TrieNode> trie_;
    std::vector<std::string> patterns_;
    int node_count_;
    
    /**
     * Add pattern to trie
     */
    void addPattern(const std::string& pattern, int pattern_index);
    
    /**
     * Build failure links (BFS)
     */
    void buildFailureLinks();
    
    /**
     * Get failure link for node
     */
    int getFailureLink(int node, char c);
    
    /**
     * Follow output links to collect all matches
     */
    void collectOutputs(int node, std::vector<int>& outputs);
};

#endif // AHO_CORASICK_H

