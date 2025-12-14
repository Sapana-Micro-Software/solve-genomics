#include "dancing_links.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <set>

DancingLinks::DancingLinks(const DLXConfig& config) : config_(config), root_(nullptr) {
    root_ = new DLXNode();                                                                      // Create root header node
    stats_ = Statistics();                                                                       // Initialize statistics
}

DancingLinks::~DancingLinks() {
    clear();                                                                                      // Clean up all nodes and reset structure
    if (root_) {                                                                                  // Check if root node exists
        delete root_;                                                                             // Delete root header node
        root_ = nullptr;                                                                          // Set root pointer to null
    }
}

void DancingLinks::buildProblem(const std::string& sequence,
                               const std::vector<std::string>& patterns,
                               int pattern_length) {
    clear();                                                                                      // Clear previous problem
    
    if (sequence.empty() || patterns.empty()) {                                                  // Validate inputs
        return;                                                                                   // Return early if invalid
    }
    
    // Calculate entropy
    stats_.entropy = calculateEntropy(sequence);                                                 // Compute sequence entropy
    
    // Build sparse matrix representation
    buildSparseMatrix(sequence, patterns);                                                        // Construct sparse matrix from sequence and patterns
    
    if (config_.low_entropy_mode && stats_.entropy < 2) {                                       // Check if low entropy mode applicable
        optimizeForLowEntropy();                                                                  // Optimize structure for low-entropy sequences
    }
}

void DancingLinks::buildSparseEntropicProblem(const std::string& sequence,
                                              const std::string& pattern) {
    clear();                                                                                      // Clear any previous problem data
    
    if (sequence.empty() || pattern.empty()) {                                                   // Validate: check if inputs are non-empty
        return;                                                                                   // Return early if inputs are invalid
    }
    
    // Calculate entropy
    stats_.entropy = calculateEntropy(sequence);                                                 // Compute entropy of DNA sequence (for sparse-entropic optimization)
    
    // Find all occurrences of pattern
    std::vector<std::string> patterns = {pattern};                                              // Create pattern vector
    std::vector<int> positions;                                                                   // Store positions where pattern occurs
    
    for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {                         // Iterate through all possible positions
        if (sequence.substr(i, pattern.length()) == pattern) {                                   // Check if pattern matches at position
            positions.push_back(static_cast<int>(i));                                            // Store matching position
        }
    }
    
    // Build exact cover: each position must be covered exactly once
    // Rows: each occurrence of pattern
    // Columns: each position in sequence
    
    int num_positions = static_cast<int>(sequence.length());                                     // Total number of positions
    int num_occurrences = static_cast<int>(positions.size());                                   // Number of pattern occurrences
    
    // Create column headers for each position
    for (int i = 0; i < num_positions; ++i) {                                                   // Iterate through all positions
        DLXNode* col = createColumnHeader(i);                                                    // Create column header
        column_headers_.push_back(col);                                                          // Store column header
    }
    
    // Create rows for each pattern occurrence
    row_to_pattern_.resize(num_occurrences);                                                     // Resize pattern mapping
    row_to_position_.resize(num_occurrences);                                                    // Resize position mapping
    
    for (int row = 0; row < num_occurrences; ++row) {                                           // Iterate through pattern occurrences
        int pos = positions[row];                                                               // Get position of this occurrence
        row_to_pattern_[row] = pattern;                                                          // Store pattern for this row
        row_to_position_[row] = pos;                                                            // Store position for this row
        
        // Add nodes for each position covered by this pattern occurrence
        for (int j = 0; j < static_cast<int>(pattern.length()); ++j) {                          // Iterate through pattern length
            int col_index = pos + j;                                                             // Calculate column index
            if (col_index < num_positions) {                                                     // Check bounds
                DLXNode* node = new DLXNode();                                                   // Create new node
                node->row = row;                                                                 // Set row index
                node->column_index = col_index;                                                  // Set column index
                node->column = column_headers_[col_index];                                      // Link to column header
                addNodeToColumn(node, column_headers_[col_index]);                               // Add node to column
                all_nodes_.emplace_back(node);                                                  // Store for cleanup
            }
        }
    }
    
    stats_.num_rows = num_occurrences;                                                           // Store number of rows
    stats_.num_columns = num_positions;                                                          // Store number of columns
    stats_.num_nodes = static_cast<int>(all_nodes_.size());                                     // Store number of nodes
    stats_.density = static_cast<double>(stats_.num_nodes) / (stats_.num_rows * stats_.num_columns); // Calculate density
}

void DancingLinks::buildSparseMatrix(const std::string& sequence,
                                     const std::vector<std::string>& patterns) {
    // Build exact cover problem:
    // Rows: each pattern occurrence at each position
    // Columns: each position in sequence (must be covered exactly once)
    
    int num_positions = static_cast<int>(sequence.length());                                     // Total positions in sequence
    int row_index = 0;                                                                           // Current row index
    
    // Create column headers
    for (int i = 0; i < num_positions; ++i) {                                                   // Iterate through all positions
        DLXNode* col = createColumnHeader(i);                                                    // Create column header
        column_headers_.push_back(col);                                                          // Store column header
    }
    
    // For each pattern, find all occurrences
    for (const std::string& pattern : patterns) {                                                // Iterate through patterns
        if (pattern.empty()) {                                                                   // Skip empty patterns
            continue;                                                                             // Continue to next pattern
        }
        
        for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {                   // Iterate through all positions
            if (sequence.substr(i, pattern.length()) == pattern) {                               // Check if pattern matches
                row_to_pattern_.push_back(pattern);                                             // Store pattern for this row
                row_to_position_.push_back(static_cast<int>(i));                                // Store position for this row
                
                // Add nodes for each position covered by this pattern
                for (size_t j = 0; j < pattern.length(); ++j) {                                 // Iterate through pattern length
                    int col_index = static_cast<int>(i + j);                                     // Calculate column index
                    if (col_index < num_positions) {                                             // Check bounds
                        DLXNode* node = new DLXNode();                                           // Create new node
                        node->row = row_index;                                                   // Set row index
                        node->column_index = col_index;                                          // Set column index
                        node->column = column_headers_[col_index];                              // Link to column header
                        addNodeToColumn(node, column_headers_[col_index]);                       // Add node to column
                        all_nodes_.emplace_back(node);                                          // Store for cleanup
                    }
                }
                row_index++;                                                                     // Increment row index
            }
        }
    }
    
    stats_.num_rows = row_index;                                                                 // Store number of rows
    stats_.num_columns = num_positions;                                                         // Store number of columns
    stats_.num_nodes = static_cast<int>(all_nodes_.size());                                    // Store number of nodes
    if (stats_.num_rows > 0 && stats_.num_columns > 0) {                                        // Check for division by zero
        stats_.density = static_cast<double>(stats_.num_nodes) / (stats_.num_rows * stats_.num_columns); // Calculate matrix density
    }
}

DancingLinks::ExactCoverSolution DancingLinks::solve() {
    ExactCoverSolution solution;                                                                 // Initialize solution structure
    auto start_time = std::chrono::high_resolution_clock::now();                              // Start timing
    
    if (!root_ || column_headers_.empty()) {                                                    // Validate problem is built
        solution.found = false;                                                                  // Mark as not found
        return solution;                                                                         // Return empty solution
    }
    
    current_solution_.clear();                                                                   // Clear current solution
    all_solutions_.clear();                                                                      // Clear all solutions
    
    bool found = search(0, current_solution_);                                                  // Start recursive search
    
    solution.found = found;                                                                      // Set found flag
    if (found && !current_solution_.empty()) {                                                  // Check if solution found
        solution.selected_rows = current_solution_;                                             // Store selected rows
        solution.num_solutions = 1;                                                              // One solution found
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();                                 // End timing
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); // Calculate duration
    solution.solve_time = duration.count();                                                     // Store solve time
    
    return solution;                                                                             // Return solution
}

bool DancingLinks::hasSolution() {
    ExactCoverSolution solution = solve();                                                       // Solve problem
    return solution.found;                                                                       // Return whether solution exists
}

std::vector<DancingLinks::ExactCoverSolution> DancingLinks::getAllSolutions() {
    std::vector<ExactCoverSolution> solutions;                                                   // Initialize solutions vector
    auto start_time = std::chrono::high_resolution_clock::now();                              // Start timing
    
    if (!root_ || column_headers_.empty()) {                                                    // Validate problem is built
        return solutions;                                                                         // Return empty vector
    }
    
    current_solution_.clear();                                                                   // Clear current solution
    all_solutions_.clear();                                                                      // Clear all solutions
    
    // Modify config temporarily to find all solutions
    bool old_find_all = config_.find_all_solutions;                                             // Save old setting
    config_.find_all_solutions = true;                                                          // Enable finding all solutions
    search(0, current_solution_);                                                               // Search for all solutions
    config_.find_all_solutions = old_find_all;                                                 // Restore old setting
    
    // Convert to solution structures
    for (const auto& sol_rows : all_solutions_) {                                               // Iterate through found solutions
        ExactCoverSolution sol;                                                                 // Create solution structure
        sol.selected_rows = sol_rows;                                                           // Store selected rows
        sol.found = true;                                                                        // Mark as found
        sol.num_solutions = static_cast<int>(all_solutions_.size());                            // Store total count
        solutions.push_back(sol);                                                                // Add to solutions
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();                                 // End timing
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); // Calculate duration
    if (!solutions.empty()) {                                                                   // Check if solutions found
        solutions[0].solve_time = duration.count();                                            // Store solve time
    }
    
    return solutions;                                                                            // Return all solutions
}

void DancingLinks::clear() {
    // Clear all nodes
    all_nodes_.clear();                                                                          // Clear all node pointers
    column_headers_.clear();                                                                     // Clear column headers
    
    // Reset statistics
    stats_ = Statistics();                                                                       // Reset statistics
    
    // Clear mappings
    row_to_pattern_.clear();                                                                     // Clear pattern mapping
    row_to_position_.clear();                                                                    // Clear position mapping
    
    // Clear solutions
    current_solution_.clear();                                                                    // Clear current solution
    all_solutions_.clear();                                                                      // Clear all solutions
}

DancingLinks::Statistics DancingLinks::getStatistics() const {
    return stats_;                                                                                // Return statistics
}

DancingLinks::DLXNode* DancingLinks::createColumnHeader(int col_index) {
    DLXNode* col = new DLXNode();                                                               // Create new column header
    col->column_index = col_index;                                                              // Set column index
    col->size = 0;                                                                               // Initialize size to zero
    
    // Link to root
    col->right = root_;                                                                          // Set right pointer to root
    col->left = root_->left;                                                                     // Set left pointer
    root_->left->right = col;                                                                    // Update root's left neighbor
    root_->left = col;                                                                           // Update root's left pointer
    
    return col;                                                                                   // Return column header
}

void DancingLinks::addNodeToColumn(DLXNode* node, DLXNode* column) {
    // Add node to bottom of column
    node->down = column;                                                                         // Set node's down pointer to column
    node->up = column->up;                                                                       // Set node's up pointer
    column->up->down = node;                                                                     // Update column's up neighbor
    column->up = node;                                                                            // Update column's up pointer
    column->size++;                                                                              // Increment column size
    
    // Link horizontally (if other nodes in same row exist)
    node->left = node;                                                                            // Initialize left pointer to self
    node->right = node;                                                                          // Initialize right pointer to self
}

void DancingLinks::coverColumn(DLXNode* column) {
    column->right->left = column->left;                                                          // Unlink column from horizontal list
    column->left->right = column->right;                                                         // Update left neighbor's right pointer
    
    // Remove all rows in this column
    DLXNode* row = column->down;                                                                 // Start from first row in column
    while (row != column) {                                                                      // Iterate through all rows
        DLXNode* node = row->right;                                                              // Start from row's right neighbor
        while (node != row) {                                                                    // Iterate through row
            node->down->up = node->up;                                                           // Unlink node from column
            node->up->down = node->down;                                                          // Update up neighbor's down pointer
            node->column->size--;                                                                 // Decrement column size
            node = node->right;                                                                  // Move to next node in row
        }
        row = row->down;                                                                          // Move to next row
    }
}

void DancingLinks::uncoverColumn(DLXNode* column) {
    // Restore all rows in this column (reverse order)
    DLXNode* row = column->up;                                                                   // Start from last row in column
    while (row != column) {                                                                      // Iterate through all rows
        DLXNode* node = row->left;                                                               // Start from row's left neighbor
        while (node != row) {                                                                    // Iterate through row
            node->column->size++;                                                                 // Increment column size
            node->down->up = node;                                                               // Restore node to column
            node->up->down = node;                                                               // Update up neighbor's down pointer
            node = node->left;                                                                   // Move to previous node in row
        }
        row = row->up;                                                                            // Move to previous row
    }
    
    column->right->left = column;                                                                 // Restore column to horizontal list
    column->left->right = column;                                                                // Update left neighbor's right pointer
}

DancingLinks::DLXNode* DancingLinks::chooseColumn() {
    DLXNode* best = nullptr;                                                                     // Best column (minimum size)
    int min_size = std::numeric_limits<int>::max();                                             // Minimum size found
    
    DLXNode* col = root_->right;                                                                 // Start from first column
    while (col != root_) {                                                                       // Iterate through all columns
        if (col->size < min_size) {                                                             // Check if smaller size
            min_size = col->size;                                                                // Update minimum size
            best = col;                                                                          // Update best column
        }
        col = col->right;                                                                        // Move to next column
    }
    
    return best;                                                                                 // Return column with minimum size
}

bool DancingLinks::search(int depth, std::vector<int>& solution) {
    if (root_->right == root_) {                                                                 // Check if all columns covered
        if (config_.find_all_solutions) {                                                       // Check if finding all solutions
            all_solutions_.push_back(solution);                                                 // Store solution
            return all_solutions_.size() >= static_cast<size_t>(config_.max_solutions);         // Check if reached max
        }
        return true;                                                                              // Solution found
    }
    
    DLXNode* col = chooseColumn();                                                               // Choose column with minimum size
    if (col == nullptr || col->size == 0) {                                                     // Check if no valid column
        return false;                                                                             // No solution
    }
    
    coverColumn(col);                                                                             // Cover chosen column
    
    // Try each row in this column
    DLXNode* row = col->down;                                                                    // Start from first row
    while (row != col) {                                                                         // Iterate through all rows
        solution.push_back(row->row);                                                            // Add row to solution
        
        // Cover all columns in this row
        DLXNode* node = row->right;                                                              // Start from row's right neighbor
        while (node != row) {                                                                    // Iterate through row
            coverColumn(node->column);                                                           // Cover column
            node = node->right;                                                                  // Move to next node
        }
        
        // Recursively search
        if (search(depth + 1, solution)) {                                                      // Recursive call
            if (!config_.find_all_solutions) {                                                  // Check if finding all solutions
                uncoverColumn(col);                                                              // Uncover column before returning
                return true;                                                                      // Solution found
            }
        }
        
        // Backtrack: remove row from solution
        solution.pop_back();                                                                     // Remove row from solution
        
        // Uncover all columns in this row
        node = row->left;                                                                        // Start from row's left neighbor
        while (node != row) {                                                                    // Iterate through row
            uncoverColumn(node->column);                                                         // Uncover column
            node = node->left;                                                                   // Move to previous node
        }
        
        row = row->down;                                                                          // Move to next row
    }
    
    uncoverColumn(col);                                                                          // Uncover column before returning
    return false;                                                                                // No solution found
}

int DancingLinks::calculateEntropy(const std::string& sequence) {
    if (sequence.empty()) {                                                                      // Check if sequence is empty
        return 0;                                                                                 // Return zero entropy
    }
    
    std::map<char, int> counts;                                                                  // Count occurrences of each character
    for (char c : sequence) {                                                                    // Iterate through sequence
        counts[c]++;                                                                              // Increment count
    }
    
    double entropy = 0.0;                                                                        // Initialize entropy
    int length = static_cast<int>(sequence.length());                                           // Sequence length
    
    for (const auto& pair : counts) {                                                            // Iterate through character counts
        double p = static_cast<double>(pair.second) / length;                                   // Calculate probability
        if (p > 0.0) {                                                                          // Check if probability is positive
            entropy -= p * std::log2(p);                                                         // Add to entropy
        }
    }
    
    return static_cast<int>(entropy * 10);                                                       // Return entropy scaled by 10
}

void DancingLinks::optimizeForLowEntropy() {
    // For low-entropy sequences, we can optimize by:
    // 1. Merging identical columns
    // 2. Reordering columns by frequency
    // 3. Using more efficient data structures
    
    // This is a placeholder for optimization
    // In a full implementation, would reorder columns and merge duplicates
}

