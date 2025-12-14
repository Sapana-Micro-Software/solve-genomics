#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>

/**
 * Suffix Tree for efficient substring search
 * Time: O(m) for construction, O(m) for search where m is pattern length
 * Space: O(n) where n is text length
 */
class SuffixTree {
public:
    SuffixTree();
    
    /**
     * Build suffix tree from text
     * @param text Text to build tree from
     */
    void build(const std::string& text);
    
    /**
     * Search for pattern in suffix tree
     * @param pattern Pattern to search for
     * @return SearchResult with all occurrences
     */
    SearchResult search(const std::string& pattern);
    
    /**
     * Check if pattern exists
     * @param pattern Pattern to check
     * @return True if pattern found
     */
    bool contains(const std::string& pattern);
    
    /**
     * Find longest common substring
     * @param other Another suffix tree
     * @return Longest common substring
     */
    std::string longestCommonSubstring(SuffixTree& other);
    
private:
    struct Node {
        int start;                    // Start position in text
        int end;                      // End position (or -1 for leaf)
        std::map<char, int> children; // Character -> child node index
        int suffix_link;              // Suffix link for Ukkonen's algorithm
        int leaf_count;               // Number of leaves in subtree
        
        Node(int s, int e) : start(s), end(e), suffix_link(-1), leaf_count(0) {}
    };
    
    std::string text_;
    std::vector<Node> nodes_;
    int root_;
    int active_node_;
    int active_edge_;
    int active_length_;
    int remaining_;
    int leaf_end_;
    
    /**
     * Create new node
     */
    int newNode(int start, int end);
    
    /**
     * Get edge length
     */
    int edgeLength(int node);
    
    /**
     * Walk down tree
     */
    bool walkDown(int node);
    
    /**
     * Extend tree (Ukkonen's algorithm)
     */
    void extend(int pos);
    
    /**
     * Build suffix tree using Ukkonen's algorithm
     */
    void buildUkkonen();
    
    /**
     * Search pattern in tree
     */
    bool searchPattern(const std::string& pattern, int& node, int& edge, int& length);
    
    /**
     * Collect all leaf positions
     */
    void collectLeaves(int node, std::vector<int>& positions);
};

#endif // SUFFIX_TREE_H

