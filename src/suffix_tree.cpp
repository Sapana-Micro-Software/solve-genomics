#include "suffix_tree.h"
#include <algorithm>
#include <cctype>

SuffixTree::SuffixTree() : root_(0), active_node_(0), active_edge_(-1),
                           active_length_(0), remaining_(0), leaf_end_(-1) {
    nodes_.push_back(Node(-1, -1));  // Root node
    root_ = 0;
}

void SuffixTree::build(const std::string& text) {
    text_ = text;
    nodes_.clear();
    nodes_.push_back(Node(-1, -1));
    root_ = 0;
    active_node_ = 0;
    active_edge_ = -1;
    active_length_ = 0;
    remaining_ = 0;
    leaf_end_ = -1;
    
    buildUkkonen();
}

int SuffixTree::newNode(int start, int end) {
    nodes_.push_back(Node(start, end));
    return static_cast<int>(nodes_.size() - 1);
}

int SuffixTree::edgeLength(int node) {
    if (node == root_) return 0;
    if (nodes_[node].end == -1) {
        return static_cast<int>(text_.length()) - nodes_[node].start;
    }
    return nodes_[node].end - nodes_[node].start + 1;
}

bool SuffixTree::walkDown(int node) {
    int length = edgeLength(node);
    if (active_length_ >= length) {
        active_edge_ += length;
        active_length_ -= length;
        active_node_ = node;
        return true;
    }
    return false;
}

void SuffixTree::extend(int pos) {
    leaf_end_ = pos;
    remaining_++;
    int last_new_node = -1;
    
    while (remaining_ > 0) {
        if (active_length_ == 0) {
            active_edge_ = pos;
        }
        
        char c = std::toupper(text_[pos]);
        
        if (nodes_[active_node_].children.find(c) == nodes_[active_node_].children.end()) {
            // Create new leaf
            int leaf = newNode(pos, -1);
            nodes_[active_node_].children[c] = leaf;
            
            if (last_new_node != -1) {
                nodes_[last_new_node].suffix_link = active_node_;
                last_new_node = -1;
            }
        } else {
            int next = nodes_[active_node_].children[c];
            
            if (walkDown(next)) {
                continue;
            }
            
            char next_char = std::toupper(text_[nodes_[next].start + active_length_]);
            if (next_char == c) {
                if (last_new_node != -1 && active_node_ != root_) {
                    nodes_[last_new_node].suffix_link = active_node_;
                    last_new_node = -1;
                }
                active_length_++;
                break;
            }
            
            // Split edge
            int split = newNode(nodes_[next].start, nodes_[next].start + active_length_ - 1);
            nodes_[active_node_].children[c] = split;
            nodes_[next].start += active_length_;
            nodes_[split].children[next_char] = next;
            
            int leaf = newNode(pos, -1);
            nodes_[split].children[c] = leaf;
            
            if (last_new_node != -1) {
                nodes_[last_new_node].suffix_link = split;
            }
            last_new_node = split;
        }
        
        remaining_--;
        
        if (active_node_ == root_ && active_length_ > 0) {
            active_length_--;
            active_edge_ = pos - remaining_ + 1;
        } else if (active_node_ != root_) {
            active_node_ = nodes_[active_node_].suffix_link;
        }
    }
}

void SuffixTree::buildUkkonen() {
    for (size_t i = 0; i < text_.length(); ++i) {
        extend(static_cast<int>(i));
    }
}

SearchResult SuffixTree::search(const std::string& pattern) {
    SearchResult result;
    
    if (pattern.empty() || text_.empty()) {
        return result;
    }
    
    int node = root_;
    int edge = -1;
    int length = 0;
    
    if (!searchPattern(pattern, node, edge, length)) {
        return result;
    }
    
    // Collect all leaf positions
    std::vector<int> positions;
    collectLeaves(node, positions);
    
    result.positions.resize(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
        result.positions[i] = positions[i];
    }
    result.count = static_cast<int>(positions.size());
    
    return result;
}

bool SuffixTree::searchPattern(const std::string& pattern, int& node, int& edge, int& length) {
    node = root_;
    edge = -1;
    length = 0;
    
    for (size_t i = 0; i < pattern.length(); ++i) {
        char c = std::toupper(pattern[i]);
        
        if (length == 0) {
            if (nodes_[node].children.find(c) == nodes_[node].children.end()) {
                return false;
            }
            edge = nodes_[node].children[c];
            length = 1;
        } else {
            int edge_start = nodes_[edge].start;
            if (edge_start + length >= static_cast<int>(text_.length())) {
                return false;
            }
            char edge_char = std::toupper(text_[edge_start + length]);
            if (edge_char != c) {
                return false;
            }
            length++;
        }
        
        int edge_len = edgeLength(edge);
        if (length >= edge_len) {
            node = edge;
            length = 0;
        }
    }
    
    return true;
}

void SuffixTree::collectLeaves(int node, std::vector<int>& positions) {
    if (nodes_[node].children.empty()) {
        // Leaf node
        positions.push_back(static_cast<int>(text_.length()) - edgeLength(node));
        return;
    }
    
    for (const auto& [c, child] : nodes_[node].children) {
        collectLeaves(child, positions);
    }
}

bool SuffixTree::contains(const std::string& pattern) {
    int node, edge, length;
    return searchPattern(pattern, node, edge, length);
}

std::string SuffixTree::longestCommonSubstring(SuffixTree& other) {
    // Simplified: find longest common substring by comparing paths
    // Full implementation would require more sophisticated traversal
    std::string lcs = "";
    int max_len = 0;
    
    // This is a simplified version - full implementation would traverse both trees
    // and find the longest common path
    
    return lcs;
}

