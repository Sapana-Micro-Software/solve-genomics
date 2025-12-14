#ifndef SKIP_GRAPH_H
#define SKIP_GRAPH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>
#include <set>
#include <functional>

/**
 * Skip-Graph for Hierarchical DNA Sequence Indexing
 * Provides logarithmic-time hash-based lookup for pre-cached patterns and subsequences
 */
class SkipGraph {
public:
    /**
     * Node in the skip-graph
     */
    struct SkipNode {
        std::string key;                                                                      // Hash key (subsequence hash)
        std::vector<int> positions;                                                            // Positions where this subsequence occurs
        std::vector<std::shared_ptr<SkipNode>> next;                                           // Next nodes at each level
        std::vector<std::shared_ptr<SkipNode>> prev;                                           // Previous nodes at each level
        int level;                                                                             // Level in the skip-graph
        
        SkipNode(const std::string& k, int lvl) : key(k), level(lvl) {                        // Constructor
            next.resize(lvl + 1, nullptr);                                                     // Initialize next pointers
            prev.resize(lvl + 1, nullptr);                                                     // Initialize prev pointers
        }
    };
    
    /**
     * Search result from skip-graph lookup
     */
    struct SkipGraphResult {
        std::vector<int> positions;                                                            // All positions found
        int num_matches;                                                                       // Number of matches
        int levels_searched;                                                                   // Number of levels searched
        double search_time;                                                                    // Search time in microseconds
        bool found;                                                                             // Whether pattern was found
        
        SkipGraphResult() : num_matches(0), levels_searched(0), search_time(0.0), found(false) {}
    };
    
    /**
     * Configuration for skip-graph construction
     */
    struct GraphConfig {
        int max_levels;                                                                        // Maximum number of levels
        int subsequence_length;                                                                 // Length of subsequences to index
        bool use_hierarchical;                                                                 // Use hierarchical indexing
        int hash_function;                                                                     // Hash function type (0=simple, 1=rolling)
        double probability;                                                                   // Probability for level promotion (0.5 for skip list)
        
        GraphConfig() : max_levels(16), subsequence_length(4), use_hierarchical(true), 
                       hash_function(1), probability(0.5) {}
    };
    
    SkipGraph(const GraphConfig& config = GraphConfig());
    ~SkipGraph();
    
    /**
     * Build skip-graph index from DNA sequence
     * @param sequence DNA sequence to index
     */
    void buildIndex(const std::string& sequence);
    
    /**
     * Pre-cache and index all subsequences of specified length
     * @param sequence DNA sequence
     * @param subsequence_length Length of subsequences to cache
     */
    void preCacheSubsequences(const std::string& sequence, int subsequence_length = -1);
    
    /**
     * Add a pattern to the index
     * @param pattern Pattern to index
     * @param position Position where pattern occurs
     */
    void addPattern(const std::string& pattern, int position);
    
    /**
     * Search for pattern using hash-based lookup
     * @param pattern Pattern to find
     * @return SkipGraphResult with positions and metadata
     */
    SkipGraphResult search(const std::string& pattern);
    
    /**
     * Search for pattern at specific level
     * @param pattern Pattern to find
     * @param level Level to search at
     * @return SkipGraphResult
     */
    SkipGraphResult searchAtLevel(const std::string& pattern, int level);
    
    /**
     * Get all positions for a hash key
     * @param hash_key Hash key to lookup
     * @return Vector of positions
     */
    std::vector<int> getPositions(const std::string& hash_key);
    
    /**
     * Check if pattern is indexed
     * @param pattern Pattern to check
     * @return True if indexed
     */
    bool isIndexed(const std::string& pattern);
    
    /**
     * Get statistics about the skip-graph
     */
    struct Statistics {
        int total_nodes;                                                                      // Total number of nodes
        int total_subsequences;                                                                // Total indexed subsequences
        int max_level;                                                                         // Maximum level used
        std::map<int, int> nodes_per_level;                                                    // Nodes at each level
        size_t memory_usage;                                                                    // Estimated memory usage
        
        Statistics() : total_nodes(0), total_subsequences(0), max_level(0), memory_usage(0) {}
    };
    
    Statistics getStatistics() const;
    
    /**
     * Clear all indexed data
     */
    void clear();
    
private:
    /**
     * Hash function for subsequences
     */
    std::string hashSubsequence(const std::string& subsequence);
    
    /**
     * Rolling hash function
     */
    std::string rollingHash(const std::string& subsequence);
    
    /**
     * Simple hash function
     */
    std::string simpleHash(const std::string& subsequence);
    
    /**
     * Generate random level for new node
     */
    int randomLevel();
    
    /**
     * Find node at specific level
     */
    std::shared_ptr<SkipNode> findNodeAtLevel(const std::string& key, int level);
    
    /**
     * Insert node into skip-graph
     */
    void insertNode(const std::string& key, int position);
    
    /**
     * Build hierarchical levels
     */
    void buildHierarchicalLevels();
    
    /**
     * Get or create node for key
     */
    std::shared_ptr<SkipNode> getOrCreateNode(const std::string& key);
    
    /**
     * Head nodes for each level
     */
    std::vector<std::shared_ptr<SkipNode>> heads_;
    
    /**
     * Hash table for quick lookup
     */
    std::unordered_map<std::string, std::shared_ptr<SkipNode>> hash_table_;
    
    /**
     * Configuration
     */
    GraphConfig config_;
    
    /**
     * Original sequence
     */
    std::string sequence_;
    
    /**
     * Indexed subsequences
     */
    std::set<std::string> indexed_subsequences_;
    
    /**
     * Random number generator
     */
    std::function<int()> random_level_generator_;
    
    /**
     * Statistics
     */
    mutable Statistics stats_;
};

#endif // SKIP_GRAPH_H

