#ifndef DANCING_LINKS_H
#define DANCING_LINKS_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <memory>
#include <functional>

/**
 * Dancing Links (Algorithm X) for Exact Cover on Sparse-Entropic DNA Vectors
 * Efficiently solves exact cover problems in low-entropy (highly repetitive) DNA sequences
 */
class DancingLinks {
public:
    /**
     * Node in the dancing links structure
     */
    struct DLXNode {
        DLXNode* left;                                                                          // Left neighbor
        DLXNode* right;                                                                          // Right neighbor
        DLXNode* up;                                                                             // Up neighbor
        DLXNode* down;                                                                           // Down neighbor
        DLXNode* column;                                                                        // Column header
        int row;                                                                                 // Row index
        int column_index;                                                                        // Column index
        int size;                                                                                // Size of column (for column headers)
        
        DLXNode() : left(this), right(this), up(this), down(this),                              // Initialize circular links
                   column(this), row(-1), column_index(-1), size(0) {}
    };
    
    /**
     * Solution to exact cover problem
     */
    struct ExactCoverSolution {
        std::vector<int> selected_rows;                                                          // Rows selected in solution
        bool found;                                                                              // Whether solution exists
        int num_solutions;                                                                      // Number of solutions found
        double solve_time;                                                                      // Time to solve (microseconds)
        
        ExactCoverSolution() : found(false), num_solutions(0), solve_time(0.0) {}
    };
    
    /**
     * Configuration for dancing links
     */
    struct DLXConfig {
        bool find_all_solutions;                                                                 // Find all solutions or just first
        int max_solutions;                                                                       // Maximum solutions to find
        bool sparse_mode;                                                                        // Optimize for sparse matrices
        bool low_entropy_mode;                                                                  // Optimize for low-entropy sequences
        
        DLXConfig() : find_all_solutions(false), max_solutions(1), 
                     sparse_mode(true), low_entropy_mode(true) {}
    };
    
    DancingLinks(const DLXConfig& config = DLXConfig());
    ~DancingLinks();
    
    /**
     * Build exact cover problem from DNA sequence and patterns
     * @param sequence DNA sequence
     * @param patterns Patterns to match
     * @param pattern_length Length of patterns (if 0, use pattern sizes)
     */
    void buildProblem(const std::string& sequence,
                     const std::vector<std::string>& patterns,
                     int pattern_length = 0);
    
    /**
     * Build exact cover for sparse-entropic vector
     * Optimized for low-entropy (highly repetitive) sequences
     * @param sequence DNA sequence (sparse-entropic)
     * @param pattern Pattern to find exact cover for
     */
    void buildSparseEntropicProblem(const std::string& sequence,
                                   const std::string& pattern);
    
    /**
     * Solve exact cover problem
     * @return ExactCoverSolution
     */
    ExactCoverSolution solve();
    
    /**
     * Check if exact cover exists
     * @return True if solution exists
     */
    bool hasSolution();
    
    /**
     * Get all solutions (up to max_solutions)
     * @return Vector of solutions
     */
    std::vector<ExactCoverSolution> getAllSolutions();
    
    /**
     * Clear problem and reset
     */
    void clear();
    
    /**
     * Get statistics about the problem
     */
    struct Statistics {
        int num_rows;                                                                            // Number of rows (constraints)
        int num_columns;                                                                         // Number of columns (elements)
        int num_nodes;                                                                           // Number of nodes in structure
        double density;                                                                          // Matrix density (0-1)
        int entropy;                                                                             // Estimated entropy of sequence
        
        Statistics() : num_rows(0), num_columns(0), num_nodes(0), density(0.0), entropy(0) {}
    };
    
    Statistics getStatistics() const;
    
private:
    /**
     * Create column header
     */
    DLXNode* createColumnHeader(int col_index);
    
    /**
     * Add node to column
     */
    void addNodeToColumn(DLXNode* node, DLXNode* column);
    
    /**
     * Remove column from structure
     */
    void coverColumn(DLXNode* column);
    
    /**
     * Restore column to structure
     */
    void uncoverColumn(DLXNode* column);
    
    /**
     * Choose column with minimum size (heuristic)
     */
    DLXNode* chooseColumn();
    
    /**
     * Recursive search for exact cover
     */
    bool search(int depth, std::vector<int>& solution);
    
    /**
     * Calculate entropy of sequence
     */
    int calculateEntropy(const std::string& sequence);
    
    /**
     * Build sparse matrix representation
     */
    void buildSparseMatrix(const std::string& sequence,
                          const std::vector<std::string>& patterns);
    
    /**
     * Optimize for low-entropy sequences
     */
    void optimizeForLowEntropy();
    
    /**
     * Root node (header)
     */
    DLXNode* root_;
    
    /**
     * Column headers
     */
    std::vector<DLXNode*> column_headers_;
    
    /**
     * All nodes (for cleanup)
     */
    std::vector<std::unique_ptr<DLXNode>> all_nodes_;
    
    /**
     * Configuration
     */
    DLXConfig config_;
    
    /**
     * Current solution being built
     */
    std::vector<int> current_solution_;
    
    /**
     * All found solutions
     */
    std::vector<std::vector<int>> all_solutions_;
    
    /**
     * Statistics
     */
    Statistics stats_;
    
    /**
     * Row to pattern mapping
     */
    std::vector<std::string> row_to_pattern_;
    
    /**
     * Row to position mapping
     */
    std::vector<int> row_to_position_;
};

#endif // DANCING_LINKS_H

