#ifndef MCMC_SEARCH_H
#define MCMC_SEARCH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <random>

/**
 * MCMC-based pattern matching for DNA sequences
 * Uses Markov Chain Monte Carlo to mutate patterns and find matches
 */
class MCMCSearch {
public:
    struct MCMCResult {
        std::vector<size_t> positions;      // Positions where matches found
        std::string evolved_pattern;       // Final evolved pattern
        int iterations;                     // Number of MCMC iterations
        double final_fitness;              // Fitness of final pattern
        int count;                         // Number of matches found
        
        MCMCResult() : iterations(0), final_fitness(0.0), count(0) {}
    };
    
    MCMCSearch(int max_iterations = 1000, 
               double mutation_rate = 0.1,
               double temperature = 1.0);
    
    /**
     * Search using MCMC to evolve pattern toward matches
     * @param sequence DNA sequence to search in
     * @param initial_pattern Starting pattern (may not match exactly)
     * @return MCMCResult with evolved pattern and matches
     */
    MCMCResult search(const std::string& sequence, const std::string& initial_pattern);
    
    /**
     * Calculate fitness of pattern (number of matches in sequence)
     */
    int calculateFitness(const std::string& sequence, const std::string& pattern);
    
    /**
     * Mutate pattern (DNA-specific mutations)
     * @param pattern Pattern to mutate
     * @param mutation_rate Probability of mutation per base
     * @return Mutated pattern
     */
    std::string mutatePattern(const std::string& pattern, double mutation_rate);
    
    /**
     * Metropolis-Hastings acceptance criterion
     * @param current_fitness Current pattern fitness
     * @param proposed_fitness Proposed pattern fitness
     * @param temperature Temperature parameter
     * @return True if proposal accepted
     */
    bool acceptProposal(int current_fitness, int proposed_fitness, double temperature);
    
    /**
     * Simulated annealing search (temperature decreases over time)
     */
    MCMCResult simulatedAnnealingSearch(const std::string& sequence, 
                                        const std::string& initial_pattern);
    
private:
    int max_iterations_;
    double mutation_rate_;
    double temperature_;
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_dist_;
    
    /**
     * DNA-specific mutation types
     */
    enum MutationType {
        SUBSTITUTION,  // A->T, G->C, etc.
        INSERTION,     // Insert random base
        DELETION       // Delete base
    };
    
    /**
     * Apply mutation to pattern
     */
    std::string applyMutation(const std::string& pattern, MutationType type, size_t position);
    
    /**
     * Get random DNA base
     */
    char randomBase();
    
    /**
     * Cool temperature (for simulated annealing)
     */
    double coolTemperature(int iteration, int max_iter);
};

#endif // MCMC_SEARCH_H

