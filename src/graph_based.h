#ifndef GRAPH_BASED_H
#define GRAPH_BASED_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>
#include <set>

/**
 * Graph-based approaches for sequence alignment and search
 */
class GraphBased {
public:
    /**
     * De Bruijn Graph for sequence assembly and pattern matching
     */
    class DeBruijnGraph {
    public:
        DeBruijnGraph(int k);
        
        /**
         * Build graph from sequences
         * @param sequences Sequences to build graph from
         */
        void build(const std::vector<std::string>& sequences);
        
        /**
         * Find Eulerian path (sequence reconstruction)
         * @return Reconstructed sequence
         */
        std::string findEulerianPath();
        
        /**
         * Search for pattern in graph
         * @param pattern Pattern to find
         * @return True if pattern found
         */
        bool search(const std::string& pattern);
        
    private:
        int k_;
        std::map<std::string, std::vector<std::string>> graph_;  // k-mer -> next k-mers
        std::map<std::string, int> in_degree_;
        std::map<std::string, int> out_degree_;
        
        /**
         * Generate k-mers from sequence
         */
        std::vector<std::string> generateKmers(const std::string& sequence);
        
        /**
         * Find starting node for Eulerian path
         */
        std::string findStartNode();
    };
    
    /**
     * Sequence Graph for alignment
     */
    class SequenceGraph {
    public:
        /**
         * Build graph from multiple sequences
         * @param sequences Sequences to align
         */
        void build(const std::vector<std::string>& sequences);
        
        /**
         * Find consensus sequence
         * @return Consensus sequence
         */
        std::string findConsensus();
        
        /**
         * Align sequences using graph
         * @return Aligned sequences
         */
        std::vector<std::string> align();
        
    private:
        struct Node {
            char base;
            std::vector<int> next_nodes;
            int count;
            
            Node(char b) : base(b), count(1) {}
        };
        
        std::vector<Node> nodes_;
        int start_node_;
        
        /**
         * Add sequence to graph
         */
        void addSequence(const std::string& sequence);
    };
    
    /**
     * Overlap Graph for sequence assembly
     */
    class OverlapGraph {
    public:
        /**
         * Build overlap graph
         * @param sequences Sequences to build graph from
         * @param min_overlap Minimum overlap required
         */
        void build(const std::vector<std::string>& sequences, int min_overlap = 3);
        
        /**
         * Find shortest superstring (approximate)
         * @return Superstring
         */
        std::string findShortestSuperstring();
        
        /**
         * Get overlap between two sequences
         * @param seq1 First sequence
         * @param seq2 Second sequence
         * @return Overlap length
         */
        int getOverlap(const std::string& seq1, const std::string& seq2);
        
    private:
        struct Edge {
            int from, to;
            int weight;  // Overlap length
            
            Edge(int f, int t, int w) : from(f), to(t), weight(w) {}
        };
        
        std::vector<std::string> sequences_;
        std::vector<std::vector<Edge>> graph_;
        
        /**
         * Greedy algorithm for shortest superstring
         */
        std::string greedySuperstring();
    };
};

#endif // GRAPH_BASED_H

