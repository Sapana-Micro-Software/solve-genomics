#include "skip_graph.h"
#include <algorithm>
#include <random>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <cstring>

SkipGraph::SkipGraph(const GraphConfig& config) : config_(config) {
    heads_.resize(config_.max_levels + 1);                                                      // Initialize head nodes for each level
    for (int i = 0; i <= config_.max_levels; ++i) {                                            // Iterate through all levels
        heads_[i] = std::make_shared<SkipNode>("", i);                                         // Create head node with empty key
    }
    
    // Initialize random level generator
    std::random_device rd;                                                                     // Random device for seeding
    std::mt19937 gen(rd());                                                                     // Mersenne Twister generator
    std::uniform_real_distribution<> dis(0.0, 1.0);                                            // Uniform distribution [0,1)
    random_level_generator_ = [this, gen, dis]() mutable -> int {                              // Lambda for random level generation
        int level = 0;                                                                          // Start at level 0
        while (level < config_.max_levels && dis(gen) < config_.probability) {                 // Promote with probability
            level++;                                                                             // Increment level
        }
        return level;                                                                            // Return generated level
    };
}

SkipGraph::~SkipGraph() {
    clear();                                                                                     // Clean up all nodes
}

void SkipGraph::buildIndex(const std::string& sequence) {
    sequence_ = sequence;                                                                       // Store original DNA sequence for reference
    clear();                                                                                     // Clear any existing index data
    
    if (sequence.empty()) {                                                                     // Validate: check if sequence is empty
        return;                                                                                  // Return early if sequence is empty
    }
    
    int subseq_len = config_.subsequence_length;                                                // Get configured subsequence length from config
    if (subseq_len <= 0) {                                                                      // Validate subsequence length
        subseq_len = 4;                                                                         // Default to length 4 if invalid
    }
    
    preCacheSubsequences(sequence, subseq_len);                                                // Pre-cache and index all subsequences of specified length
}

void SkipGraph::preCacheSubsequences(const std::string& sequence, int subsequence_length) {
    if (sequence.empty()) {                                                                     // Validate: check if sequence is empty
        return;                                                                                  // Return early if sequence is empty
    }
    
    int subseq_len = subsequence_length > 0 ? subsequence_length : config_.subsequence_length; // Use provided length or fall back to configured length
    
    for (int i = 0; i <= static_cast<int>(sequence.length()) - subseq_len; ++i) {              // Iterate through all valid starting positions for subsequences
        std::string subsequence = sequence.substr(i, subseq_len);                               // Extract subsequence of specified length
        addPattern(subsequence, i);                                                              // Add pattern to skip-graph index with its position
    }
    
    if (config_.use_hierarchical) {                                                            // Check if hierarchical indexing is enabled in configuration
        buildHierarchicalLevels();                                                              // Build hierarchical skip-graph structure with multiple levels
    }
}

void SkipGraph::addPattern(const std::string& pattern, int position) {
    if (pattern.empty()) {                                                                      // Check if pattern is empty
        return;                                                                                  // Return early if empty
    }
    
    std::string hash_key = hashSubsequence(pattern);                                            // Compute hash key for pattern
    indexed_subsequences_.insert(pattern);                                                      // Track indexed subsequence
    
    insertNode(hash_key, position);                                                             // Insert node into skip-graph
}

SkipGraph::SkipGraphResult SkipGraph::search(const std::string& pattern) {
    SkipGraphResult result;                                                                     // Initialize result structure for search
    auto start_time = std::chrono::high_resolution_clock::now();                              // Record start time for performance measurement
    
    if (pattern.empty()) {                                                                      // Validate: check if pattern is empty
        result.found = false;                                                                   // Mark as not found
        return result;                                                                           // Return empty result immediately
    }
    
    std::string hash_key = hashSubsequence(pattern);                                            // Compute hash key for pattern using configured hash function
    
    // Try hash table lookup first (O(1) average case)
    auto it = hash_table_.find(hash_key);                                                       // Perform O(1) hash table lookup
    if (it != hash_table_.end()) {                                                              // Check if hash key found in hash table
        result.positions = it->second->positions;                                             // Retrieve all positions where pattern occurs
        result.num_matches = static_cast<int>(result.positions.size());                         // Count total number of matches
        result.found = result.num_matches > 0;                                                 // Mark as found if at least one match exists
        result.levels_searched = 1;                                                             // Only searched hash table (single level)
    } else {
        // Fall back to skip-graph search (O(log n) worst case)
        int max_level = static_cast<int>(heads_.size()) - 1;                                   // Get maximum level in skip-graph
        for (int level = max_level; level >= 0; --level) {                                     // Search from top level down to base level
            std::shared_ptr<SkipNode> node = findNodeAtLevel(hash_key, level);                  // Find node at this level using key comparison
            if (node && node->key == hash_key) {                                                // Check if exact match found (key matches)
                result.positions = node->positions;                                            // Retrieve all positions from matching node
                result.num_matches = static_cast<int>(result.positions.size());                 // Count total number of matches
                result.found = result.num_matches > 0;                                         // Mark as found if matches exist
                result.levels_searched = max_level - level + 1;                                // Count number of levels searched
                break;                                                                          // Exit loop early when match found
            }
            result.levels_searched++;                                                            // Increment levels searched counter
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();                                 // Record end time
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); // Calculate search duration
    result.search_time = duration.count();                                                      // Store search time in microseconds
    
    return result;                                                                               // Return search result with positions and metadata
}

SkipGraph::SkipGraphResult SkipGraph::searchAtLevel(const std::string& pattern, int level) {
    SkipGraphResult result;                                                                     // Initialize result structure
    auto start_time = std::chrono::high_resolution_clock::now();                              // Start timing
    
    if (pattern.empty() || level < 0 || level >= static_cast<int>(heads_.size())) {            // Validate inputs
        result.found = false;                                                                   // Mark as not found
        return result;                                                                           // Return empty result
    }
    
    std::string hash_key = hashSubsequence(pattern);                                            // Compute hash key for pattern
    std::shared_ptr<SkipNode> node = findNodeAtLevel(hash_key, level);                          // Find node at specified level
    
    if (node && node->key == hash_key) {                                                        // Check if exact match found
        result.positions = node->positions;                                                    // Get positions
        result.num_matches = static_cast<int>(result.positions.size());                         // Count matches
        result.found = result.num_matches > 0;                                                 // Mark as found
    } else {
        result.found = false;                                                                   // Mark as not found
    }
    
    result.levels_searched = 1;                                                                 // Searched only one level
    auto end_time = std::chrono::high_resolution_clock::now();                                 // End timing
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); // Calculate duration
    result.search_time = duration.count();                                                      // Store search time
    
    return result;                                                                               // Return search result
}

std::vector<int> SkipGraph::getPositions(const std::string& hash_key) {
    auto it = hash_table_.find(hash_key);                                                       // Lookup in hash table
    if (it != hash_table_.end()) {                                                              // Check if found
        return it->second->positions;                                                           // Return positions
    }
    return std::vector<int>();                                                                   // Return empty vector if not found
}

bool SkipGraph::isIndexed(const std::string& pattern) {
    return indexed_subsequences_.find(pattern) != indexed_subsequences_.end();                 // Check if pattern is in indexed set
}

SkipGraph::Statistics SkipGraph::getStatistics() const {
    Statistics stats;                                                                            // Initialize statistics
    
    // Count nodes
    stats.total_nodes = static_cast<int>(hash_table_.size());                                  // Count nodes in hash table
    stats.total_subsequences = static_cast<int>(indexed_subsequences_.size());                 // Count indexed subsequences
    
    // Find maximum level
    int max_lvl = 0;                                                                            // Initialize max level
    for (const auto& pair : hash_table_) {                                                      // Iterate through all nodes
        if (pair.second) {                                                                      // Check if node exists
            max_lvl = std::max(max_lvl, pair.second->level);                                   // Update max level
            stats.nodes_per_level[pair.second->level]++;                                      // Count nodes per level
        }
    }
    stats.max_level = max_lvl;                                                                  // Store maximum level
    
    // Estimate memory usage
    stats.memory_usage = hash_table_.size() * sizeof(SkipNode) +                               // Node memory
                        indexed_subsequences_.size() * 32;                                     // Subsequence storage (approx)
    
    return stats;                                                                                // Return statistics
}

void SkipGraph::clear() {
    hash_table_.clear();                                                                        // Clear hash table
    indexed_subsequences_.clear();                                                              // Clear indexed subsequences
    
    // Reset head nodes
    for (int i = 0; i < static_cast<int>(heads_.size()); ++i) {                                 // Iterate through all levels
        heads_[i] = std::make_shared<SkipNode>("", i);                                         // Recreate head nodes
    }
    
    sequence_.clear();                                                                          // Clear stored sequence
}

std::string SkipGraph::hashSubsequence(const std::string& subsequence) {
    if (config_.hash_function == 0) {                                                          // Check hash function type
        return simpleHash(subsequence);                                                         // Use simple hash
    } else {
        return rollingHash(subsequence);                                                        // Use rolling hash
    }
}

std::string SkipGraph::rollingHash(const std::string& subsequence) {
    const int base = 256;                                                                       // Base for polynomial rolling hash (256 for ASCII)
    const int mod = 1000000007;                                                                 // Large prime modulus to prevent overflow
    
    long long hash = 0;                                                                         // Initialize hash value to zero
    long long power = 1;                                                                        // Initialize power of base (base^0 = 1)
    
    for (size_t i = 0; i < subsequence.length(); ++i) {                                         // Iterate through each character in subsequence
        hash = (hash + (subsequence[i] * power) % mod) % mod;                                  // Add character's contribution weighted by position
        power = (power * base) % mod;                                                           // Update power: base^i for next iteration
    }
    
    std::ostringstream oss;                                                                     // String stream for hexadecimal conversion
    oss << std::hex << std::setfill('0') << std::setw(16) << hash;                             // Convert hash to 16-character hex string
    return oss.str();                                                                            // Return hexadecimal hash string
}

std::string SkipGraph::simpleHash(const std::string& subsequence) {
    std::hash<std::string> hasher;                                                              // Standard library string hash function
    size_t hash_value = hasher(subsequence);                                                    // Compute hash value for subsequence
    std::ostringstream oss;                                                                     // String stream for hexadecimal conversion
    oss << std::hex << hash_value;                                                              // Convert hash value to hexadecimal string
    return oss.str();                                                                            // Return hexadecimal hash string
}

int SkipGraph::randomLevel() {
    return random_level_generator_();                                                            // Generate random level
}

std::shared_ptr<SkipGraph::SkipNode> SkipGraph::findNodeAtLevel(const std::string& key, int level) {
    if (level < 0 || level >= static_cast<int>(heads_.size())) {                                // Validate level bounds
        return nullptr;                                                                          // Return null pointer if level is invalid
    }
    
    std::shared_ptr<SkipNode> current = heads_[level];                                           // Start traversal from head node at specified level
    
    while (current->next[level] != nullptr && current->next[level]->key < key) {                // Traverse while next node's key is less than search key
        current = current->next[level];                                                          // Move to next node at this level
    }
    
    return current;                                                                              // Return node (may be exact match or predecessor node)
}

void SkipGraph::insertNode(const std::string& key, int position) {
    std::shared_ptr<SkipNode> node = getOrCreateNode(key);                                     // Get existing node or create new node for hash key
    node->positions.push_back(position);                                                        // Add position to node's position list
    
    int node_level = node->level;                                                                // Get node's maximum level
    
    // Insert node at each level from 0 to node_level
    for (int level = 0; level <= node_level; ++level) {                                         // Iterate through all levels where node exists
        std::shared_ptr<SkipNode> current = heads_[level];                                     // Start from head node at this level
        
        // Find insertion point (maintain sorted order by key)
        while (current->next[level] != nullptr && current->next[level]->key < key) {          // Traverse to find correct position
            current = current->next[level];                                                     // Move to next node
        }
        
        // Insert node into linked list at this level
        node->next[level] = current->next[level];                                               // Set node's next pointer to current's next
        if (current->next[level] != nullptr) {                                                 // Check if next node exists
            current->next[level]->prev[level] = node;                                           // Update next node's previous pointer
        }
        current->next[level] = node;                                                            // Update current's next pointer to new node
        node->prev[level] = current;                                                            // Set node's previous pointer to current
    }
}

void SkipGraph::buildHierarchicalLevels() {
    // Rebuild levels based on existing nodes
    // This creates a hierarchical structure where higher levels contain fewer nodes
    
    std::vector<std::shared_ptr<SkipNode>> all_nodes;                                          // Collect all nodes
    for (const auto& pair : hash_table_) {                                                     // Iterate through hash table
        all_nodes.push_back(pair.second);                                                       // Add node to list
    }
    
    // Sort nodes by key
    std::sort(all_nodes.begin(), all_nodes.end(),                                               // Sort nodes
              [](const std::shared_ptr<SkipNode>& a, const std::shared_ptr<SkipNode>& b) {     // Comparison function
                  return a->key < b->key;                                                      // Compare by key
              });
    
    // Rebuild skip-graph structure
    for (int level = 0; level <= config_.max_levels; ++level) {                                // Iterate through all levels
        std::shared_ptr<SkipNode> current = heads_[level];                                     // Start from head
        
        for (size_t i = 0; i < all_nodes.size(); ++i) {                                        // Iterate through sorted nodes
            if (all_nodes[i]->level >= level) {                                                // Check if node exists at this level
                current->next[level] = all_nodes[i];                                           // Link node
                all_nodes[i]->prev[level] = current;                                           // Set back link
                all_nodes[i]->next[level] = nullptr;                                           // Initialize next
                current = all_nodes[i];                                                          // Move to next node
            }
        }
    }
}

std::shared_ptr<SkipGraph::SkipNode> SkipGraph::getOrCreateNode(const std::string& key) {
    auto it = hash_table_.find(key);                                                            // Lookup node in hash table using hash key
    if (it != hash_table_.end()) {                                                              // Check if node already exists
        return it->second;                                                                       // Return existing node (reuse)
    }
    
    // Create new node for this hash key
    int level = randomLevel();                                                                  // Generate random level using probability distribution
    std::shared_ptr<SkipNode> node = std::make_shared<SkipNode>(key, level);                   // Create new skip-node with key and level
    hash_table_[key] = node;                                                                    // Store node in hash table for O(1) lookup
    
    return node;                                                                                 // Return newly created node
}

