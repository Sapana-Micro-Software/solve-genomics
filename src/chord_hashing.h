#ifndef CHORD_HASHING_H
#define CHORD_HASHING_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>

/**
 * Chord Distributed Hash Table (DHT) for DNA sequence storage and lookup
 * Based on: "Chord: A Scalable Peer-to-peer Lookup Service for Internet Applications"
 * by Ion Stoica, Robert Morris, David Karger, Frans Kaashoek, and Hari Balakrishnan
 */
class ChordHashing {
public:
    /**
     * Node in Chord ring
     */
    struct ChordNode {
        size_t node_id;                    // Node identifier (hash)
        std::string node_address;          // Node address/identifier
        std::map<size_t, std::string> finger_table;  // Finger table for routing
        std::map<size_t, std::string> data_store;    // Key-value pairs stored at this node
        
        ChordNode(size_t id, const std::string& addr) 
            : node_id(id), node_address(addr) {}
    };
    
    /**
     * Result of Chord lookup
     */
    struct LookupResult {
        size_t key_hash;                   // Hash of the key
        size_t responsible_node_id;        // Node responsible for this key
        std::string responsible_node_addr; // Address of responsible node
        std::string value;                 // Value stored at key (if found)
        int hops;                          // Number of hops to find node
        
        LookupResult() : key_hash(0), responsible_node_id(0), hops(0) {}
    };
    
    ChordHashing(int m_bits = 160);  // m-bit identifier space (default 160 for SHA-1)
    
    /**
     * Add node to Chord ring
     * @param node_address Node identifier/address
     * @return Node ID (hash)
     */
    size_t addNode(const std::string& node_address);
    
    /**
     * Remove node from Chord ring
     * @param node_id Node ID to remove
     */
    void removeNode(size_t node_id);
    
    /**
     * Hash a key (sequence) to identifier space
     * @param key Key to hash (DNA sequence)
     * @return Hash value in [0, 2^m)
     */
    size_t hashKey(const std::string& key);
    
    /**
     * Find successor node for a given key
     * @param key_hash Hash of the key
     * @return LookupResult with responsible node
     */
    LookupResult findSuccessor(size_t key_hash);
    
    /**
     * Find predecessor node for a given key
     * @param key_hash Hash of the key
     * @return Node ID of predecessor
     */
    size_t findPredecessor(size_t key_hash);
    
    /**
     * Store key-value pair in Chord DHT
     * @param key DNA sequence (key)
     * @param value Associated value/metadata
     * @return LookupResult with storage location
     */
    LookupResult store(const std::string& key, const std::string& value);
    
    /**
     * Lookup value for a key
     * @param key DNA sequence (key)
     * @return LookupResult with value
     */
    LookupResult lookup(const std::string& key);
    
    /**
     * Build finger table for a node
     * @param node_id Node to build finger table for
     */
    void buildFingerTable(size_t node_id);
    
    /**
     * Get all nodes in the ring
     * @return Vector of node IDs
     */
    std::vector<size_t> getNodes();
    
    /**
     * Get data stored at a node
     * @param node_id Node ID
     * @return Map of key hashes to values
     */
    std::map<size_t, std::string> getNodeData(size_t node_id);
    
    /**
     * Stabilize Chord ring (fix finger tables and successors)
     */
    void stabilize();
    
    /**
     * Get statistics about the Chord ring
     */
    struct RingStatistics {
        int num_nodes;
        int total_keys;
        double avg_keys_per_node;
        int max_keys_per_node;
        int min_keys_per_node;
        double avg_finger_table_size;
        
        RingStatistics() : num_nodes(0), total_keys(0), avg_keys_per_node(0.0),
                          max_keys_per_node(0), min_keys_per_node(0), 
                          avg_finger_table_size(0.0) {}
    };
    
    RingStatistics getStatistics();
    
private:
    int m_bits_;                           // Number of bits in identifier space
    size_t max_id_;                        // Maximum ID (2^m - 1)
    std::map<size_t, std::shared_ptr<ChordNode>> nodes_;  // All nodes in ring
    
    /**
     * SHA-1 like hash function (simplified for DNA sequences)
     */
    size_t sha1Hash(const std::string& input);
    
    /**
     * Check if key is in range (start, end] on Chord ring
     */
    bool inRange(size_t key, size_t start, size_t end);
    
    /**
     * Find closest preceding finger for a key
     */
    size_t closestPrecedingFinger(size_t node_id, size_t key_hash);
    
    /**
     * Get successor of a node
     */
    size_t getSuccessor(size_t node_id);
    
    /**
     * Get predecessor of a node
     */
    size_t getPredecessor(size_t node_id);
    
    /**
     * Sort nodes by ID
     */
    void sortNodes();
};

#endif // CHORD_HASHING_H

