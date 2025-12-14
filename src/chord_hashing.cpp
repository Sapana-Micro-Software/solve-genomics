#include "chord_hashing.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream>
#include <climits>

ChordHashing::ChordHashing(int m_bits) : m_bits_(m_bits) {
    if (m_bits_ < 1) m_bits_ = 160;
    if (m_bits_ > 64) m_bits_ = 64;  // Limit to 64 bits for size_t
    
    max_id_ = (1ULL << m_bits_) - 1;
}

size_t ChordHashing::sha1Hash(const std::string& input) {
    // Simplified hash function (not true SHA-1, but similar properties)
    std::hash<std::string> hasher;
    size_t hash = hasher(input);
    
    // Reduce to m_bits_ bits
    return hash & max_id_;
}

size_t ChordHashing::hashKey(const std::string& key) {
    return sha1Hash(key);
}

size_t ChordHashing::addNode(const std::string& node_address) {
    size_t node_id = hashKey(node_address);
    
    // If node already exists, return existing ID
    if (nodes_.find(node_id) != nodes_.end()) {
        return node_id;
    }
    
    // Create new node
    auto node = std::make_shared<ChordNode>(node_id, node_address);
    nodes_[node_id] = node;
    
    // Build finger table
    buildFingerTable(node_id);
    
    // Stabilize ring
    stabilize();
    
    return node_id;
}

void ChordHashing::removeNode(size_t node_id) {
    if (nodes_.find(node_id) == nodes_.end()) {
        return;
    }
    
    // Redistribute data (simplified - in practice would move to successor)
    auto node = nodes_[node_id];
    size_t successor = getSuccessor(node_id);
    
    if (successor != node_id && nodes_.find(successor) != nodes_.end()) {
        // Move data to successor
        for (const auto& [key_hash, value] : node->data_store) {
            nodes_[successor]->data_store[key_hash] = value;
        }
    }
    
    // Remove node
    nodes_.erase(node_id);
    
    // Rebuild finger tables
    stabilize();
}

bool ChordHashing::inRange(size_t key, size_t start, size_t end) {
    if (start < end) {
        return key > start && key <= end;
    } else {
        // Wraps around
        return key > start || key <= end;
    }
}

size_t ChordHashing::getSuccessor(size_t node_id) {
    if (nodes_.empty()) {
        return node_id;
    }
    
    std::vector<size_t> node_ids;
    for (const auto& [id, _] : nodes_) {
        node_ids.push_back(id);
    }
    std::sort(node_ids.begin(), node_ids.end());
    
    auto it = std::upper_bound(node_ids.begin(), node_ids.end(), node_id);
    if (it != node_ids.end()) {
        return *it;
    }
    
    // Wrap around to first node
    return node_ids.empty() ? node_id : node_ids[0];
}

size_t ChordHashing::getPredecessor(size_t node_id) {
    if (nodes_.empty()) {
        return node_id;
    }
    
    std::vector<size_t> node_ids;
    for (const auto& [id, _] : nodes_) {
        node_ids.push_back(id);
    }
    std::sort(node_ids.begin(), node_ids.end());
    
    auto it = std::lower_bound(node_ids.begin(), node_ids.end(), node_id);
    if (it != node_ids.begin()) {
        return *(--it);
    }
    
    // Wrap around to last node
    return node_ids.empty() ? node_id : node_ids.back();
}

size_t ChordHashing::findPredecessor(size_t key_hash) {
    if (nodes_.empty()) {
        return 0;
    }
    
    size_t current = nodes_.begin()->first;
    
    while (true) {
        size_t successor = getSuccessor(current);
        
        if (successor == current || inRange(key_hash, current, successor)) {
            return current;
        }
        
        // Use finger table for faster lookup
        size_t next = closestPrecedingFinger(current, key_hash);
        if (next == current) {
            current = successor;
        } else {
            current = next;
        }
    }
}

size_t ChordHashing::closestPrecedingFinger(size_t node_id, size_t key_hash) {
    if (nodes_.find(node_id) == nodes_.end()) {
        return node_id;
    }
    
    auto node = nodes_[node_id];
    
    // Check finger table in reverse order
    for (int i = m_bits_ - 1; i >= 0; --i) {
        size_t finger_id = (node_id + (1ULL << i)) & max_id_;
        
        if (node->finger_table.find(finger_id) != node->finger_table.end()) {
            size_t finger_node_id = hashKey(node->finger_table[finger_id]);
            
            if (inRange(finger_node_id, node_id, key_hash)) {
                return finger_node_id;
            }
        }
    }
    
    return node_id;
}

ChordHashing::LookupResult ChordHashing::findSuccessor(size_t key_hash) {
    LookupResult result;
    result.key_hash = key_hash;
    
    if (nodes_.empty()) {
        return result;
    }
    
    size_t predecessor = findPredecessor(key_hash);
    size_t successor = getSuccessor(predecessor);
    
    result.responsible_node_id = successor;
    if (nodes_.find(successor) != nodes_.end()) {
        result.responsible_node_addr = nodes_[successor]->node_address;
    }
    
    // Count hops (simplified)
    result.hops = 1;
    
    return result;
}

void ChordHashing::buildFingerTable(size_t node_id) {
    if (nodes_.find(node_id) == nodes_.end()) {
        return;
    }
    
    auto node = nodes_[node_id];
    node->finger_table.clear();
    
    for (int i = 0; i < m_bits_; ++i) {
        size_t finger_id = (node_id + (1ULL << i)) & max_id_;
        size_t finger_successor = findSuccessor(finger_id).responsible_node_id;
        
        if (nodes_.find(finger_successor) != nodes_.end()) {
            node->finger_table[finger_id] = nodes_[finger_successor]->node_address;
        }
    }
}

ChordHashing::LookupResult ChordHashing::store(const std::string& key, const std::string& value) {
    size_t key_hash = hashKey(key);
    LookupResult result = findSuccessor(key_hash);
    
    if (nodes_.find(result.responsible_node_id) != nodes_.end()) {
        nodes_[result.responsible_node_id]->data_store[key_hash] = value;
        result.value = value;
    }
    
    return result;
}

ChordHashing::LookupResult ChordHashing::lookup(const std::string& key) {
    size_t key_hash = hashKey(key);
    LookupResult result = findSuccessor(key_hash);
    
    if (nodes_.find(result.responsible_node_id) != nodes_.end()) {
        auto node = nodes_[result.responsible_node_id];
        if (node->data_store.find(key_hash) != node->data_store.end()) {
            result.value = node->data_store[key_hash];
        }
    }
    
    return result;
}

void ChordHashing::stabilize() {
    // Rebuild all finger tables
    for (const auto& [node_id, node] : nodes_) {
        buildFingerTable(node_id);
    }
}

std::vector<size_t> ChordHashing::getNodes() {
    std::vector<size_t> node_ids;
    for (const auto& [id, _] : nodes_) {
        node_ids.push_back(id);
    }
    std::sort(node_ids.begin(), node_ids.end());
    return node_ids;
}

std::map<size_t, std::string> ChordHashing::getNodeData(size_t node_id) {
    if (nodes_.find(node_id) != nodes_.end()) {
        return nodes_[node_id]->data_store;
    }
    return {};
}

ChordHashing::RingStatistics ChordHashing::getStatistics() {
    RingStatistics stats;
    
    if (nodes_.empty()) {
        return stats;
    }
    
    stats.num_nodes = static_cast<int>(nodes_.size());
    
    int total_keys = 0;
    int max_keys = 0;
    int min_keys = INT_MAX;
    int total_fingers = 0;
    
    for (const auto& [node_id, node] : nodes_) {
        int keys = static_cast<int>(node->data_store.size());
        total_keys += keys;
        max_keys = std::max(max_keys, keys);
        min_keys = std::min(min_keys, keys);
        total_fingers += static_cast<int>(node->finger_table.size());
    }
    
    stats.total_keys = total_keys;
    stats.avg_keys_per_node = static_cast<double>(total_keys) / stats.num_nodes;
    stats.max_keys_per_node = max_keys;
    stats.min_keys_per_node = (min_keys == INT_MAX) ? 0 : min_keys;
    stats.avg_finger_table_size = static_cast<double>(total_fingers) / stats.num_nodes;
    
    return stats;
}

void ChordHashing::sortNodes() {
    // Nodes are stored in map, already sorted by key
    // This is a placeholder for any additional sorting logic
}

