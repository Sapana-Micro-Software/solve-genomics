#include "graph_based.h"
#include <algorithm>
#include <cctype>
#include <queue>
#include <stack>

// De Bruijn Graph Implementation

GraphBased::DeBruijnGraph::DeBruijnGraph(int k) : k_(k) {
    if (k_ < 2) k_ = 3;
}

std::vector<std::string> GraphBased::DeBruijnGraph::generateKmers(const std::string& sequence) {
    std::vector<std::string> kmers;
    
    if (static_cast<int>(sequence.length()) < k_) {
        return kmers;
    }
    
    for (size_t i = 0; i <= sequence.length() - k_; ++i) {
        std::string kmer = sequence.substr(i, k_);
        for (char& c : kmer) {
            c = std::toupper(c);
        }
        kmers.push_back(kmer);
    }
    
    return kmers;
}

void GraphBased::DeBruijnGraph::build(const std::vector<std::string>& sequences) {
    graph_.clear();
    in_degree_.clear();
    out_degree_.clear();
    
    for (const std::string& sequence : sequences) {
        std::vector<std::string> kmers = generateKmers(sequence);
        
        for (size_t i = 0; i + 1 < kmers.size(); ++i) {
            std::string from = kmers[i];
            std::string to = kmers[i + 1];
            
            // Get prefix and suffix
            std::string from_suffix = from.substr(1);
            std::string to_prefix = to.substr(0, k_ - 1);
            
            if (from_suffix == to_prefix) {
                graph_[from].push_back(to);
                out_degree_[from]++;
                in_degree_[to]++;
                
                if (in_degree_.find(from) == in_degree_.end()) {
                    in_degree_[from] = 0;
                }
                if (out_degree_.find(to) == out_degree_.end()) {
                    out_degree_[to] = 0;
                }
            }
        }
    }
}

bool GraphBased::DeBruijnGraph::search(const std::string& pattern) {
    if (static_cast<int>(pattern.length()) < k_) {
        return false;
    }
    
    std::vector<std::string> pattern_kmers = generateKmers(pattern);
    
    if (pattern_kmers.empty()) {
        return false;
    }
    
    // Check if all k-mers exist in graph
    for (const std::string& kmer : pattern_kmers) {
        if (graph_.find(kmer) == graph_.end() && 
            in_degree_.find(kmer) == in_degree_.end()) {
            return false;
        }
    }
    
    // Check if k-mers form a path
    for (size_t i = 0; i + 1 < pattern_kmers.size(); ++i) {
        const std::string& from = pattern_kmers[i];
        const std::string& to = pattern_kmers[i + 1];
        
        if (graph_.find(from) == graph_.end()) {
            return false;
        }
        
        bool found = false;
        for (const std::string& next : graph_[from]) {
            if (next == to) {
                found = true;
                break;
            }
        }
        
        if (!found) {
            return false;
        }
    }
    
    return true;
}

std::string GraphBased::DeBruijnGraph::findStartNode() {
    // Find node with out_degree > in_degree
    for (const auto& [node, out] : out_degree_) {
        int in = (in_degree_.find(node) != in_degree_.end()) ? in_degree_[node] : 0;
        if (out > in) {
            return node;
        }
    }
    
    // If all balanced, return first node
    if (!graph_.empty()) {
        return graph_.begin()->first;
    }
    
    return "";
}

std::string GraphBased::DeBruijnGraph::findEulerianPath() {
    if (graph_.empty()) {
        return "";
    }
    
    std::string start = findStartNode();
    if (start.empty()) {
        return "";
    }
    
    // Hierholzer's algorithm for Eulerian path
    std::stack<std::string> stack;
    std::vector<std::string> path;
    std::map<std::string, size_t> next_edge;
    
    for (const auto& [node, _] : graph_) {
        next_edge[node] = 0;
    }
    
    stack.push(start);
    
    while (!stack.empty()) {
        std::string current = stack.top();
        
        if (next_edge[current] < graph_[current].size()) {
            std::string next = graph_[current][next_edge[current]];
            next_edge[current]++;
            stack.push(next);
        } else {
            path.push_back(current);
            stack.pop();
        }
    }
    
    std::reverse(path.begin(), path.end());
    
    // Reconstruct sequence from path
    if (path.empty()) {
        return "";
    }
    
    std::string result = path[0];
    for (size_t i = 1; i < path.size(); ++i) {
        if (!path[i].empty()) {
            result += path[i].back();
        }
    }
    
    return result;
}

// Sequence Graph Implementation

void GraphBased::SequenceGraph::build(const std::vector<std::string>& sequences) {
    nodes_.clear();
    start_node_ = -1;
    
    for (const std::string& sequence : sequences) {
        addSequence(sequence);
    }
}

void GraphBased::SequenceGraph::addSequence(const std::string& sequence) {
    if (start_node_ == -1) {
        // First sequence
        for (char c : sequence) {
            int node_id = static_cast<int>(nodes_.size());
            nodes_.push_back(Node(std::toupper(c)));
            if (node_id > 0) {
                nodes_[node_id - 1].next_nodes.push_back(node_id);
            }
        }
        start_node_ = 0;
    } else {
        // Merge with existing graph (simplified)
        // In full implementation, would use more sophisticated alignment
        size_t pos = 0;
        int current = start_node_;
        
        while (pos < sequence.length() && current < static_cast<int>(nodes_.size())) {
            char c = std::toupper(sequence[pos]);
            if (nodes_[current].base == c) {
                nodes_[current].count++;
                if (!nodes_[current].next_nodes.empty()) {
                    current = nodes_[current].next_nodes[0];
                }
                pos++;
            } else {
                break;
            }
        }
        
        // Add remaining characters
        while (pos < sequence.length()) {
            int new_node = static_cast<int>(nodes_.size());
            nodes_.push_back(Node(std::toupper(sequence[pos])));
            if (current >= 0 && current < static_cast<int>(nodes_.size()) - 1) {
                nodes_[current].next_nodes.push_back(new_node);
            }
            current = new_node;
            pos++;
        }
    }
}

std::string GraphBased::SequenceGraph::findConsensus() {
    if (nodes_.empty() || start_node_ == -1) {
        return "";
    }
    
    std::string consensus;
    int current = start_node_;
    
    while (current >= 0 && current < static_cast<int>(nodes_.size())) {
        consensus += nodes_[current].base;
        
        if (nodes_[current].next_nodes.empty()) {
            break;
        }
        
        // Choose most frequent path (simplified)
        current = nodes_[current].next_nodes[0];
    }
    
    return consensus;
}

std::vector<std::string> GraphBased::SequenceGraph::align() {
    // Simplified alignment - returns consensus
    std::string consensus = findConsensus();
    return {consensus};
}

// Overlap Graph Implementation

int GraphBased::OverlapGraph::getOverlap(const std::string& seq1, const std::string& seq2) {
    int max_overlap = 0;
    int min_len = static_cast<int>(std::min(seq1.length(), seq2.length()));
    
    for (int len = min_len; len >= 1; --len) {
        std::string suffix = seq1.substr(seq1.length() - len);
        std::string prefix = seq2.substr(0, len);
        
        bool match = true;
        for (int i = 0; i < len; ++i) {
            if (std::toupper(suffix[i]) != std::toupper(prefix[i])) {
                match = false;
                break;
            }
        }
        
        if (match) {
            max_overlap = len;
            break;
        }
    }
    
    return max_overlap;
}

void GraphBased::OverlapGraph::build(const std::vector<std::string>& sequences, int min_overlap) {
    sequences_ = sequences;
    graph_.clear();
    graph_.resize(sequences.size());
    
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (size_t j = 0; j < sequences.size(); ++j) {
            if (i != j) {
                int overlap = getOverlap(sequences[i], sequences[j]);
                if (overlap >= min_overlap) {
                    graph_[i].push_back(Edge(static_cast<int>(i), static_cast<int>(j), overlap));
                }
            }
        }
    }
}

std::string GraphBased::OverlapGraph::greedySuperstring() {
    if (sequences_.empty()) {
        return "";
    }
    
    std::vector<bool> used(sequences_.size(), false);
    std::string result = sequences_[0];
    used[0] = true;
    
    while (true) {
        int best_i = -1;
        int best_overlap = 0;
        bool append_front = false;
        
        for (size_t i = 0; i < sequences_.size(); ++i) {
            if (used[i]) continue;
            
            // Check overlap at end
            int overlap_end = getOverlap(result, sequences_[i]);
            if (overlap_end > best_overlap) {
                best_overlap = overlap_end;
                best_i = static_cast<int>(i);
                append_front = false;
            }
            
            // Check overlap at beginning
            int overlap_begin = getOverlap(sequences_[i], result);
            if (overlap_begin > best_overlap) {
                best_overlap = overlap_begin;
                best_i = static_cast<int>(i);
                append_front = true;
            }
        }
        
        if (best_i == -1) {
            break;
        }
        
        if (append_front) {
            result = sequences_[best_i] + result.substr(best_overlap);
        } else {
            result = result + sequences_[best_i].substr(best_overlap);
        }
        
        used[best_i] = true;
    }
    
    return result;
}

std::string GraphBased::OverlapGraph::findShortestSuperstring() {
    return greedySuperstring();
}

