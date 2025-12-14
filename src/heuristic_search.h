#ifndef HEURISTIC_SEARCH_H
#define HEURISTIC_SEARCH_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <queue>
#include <map>

/**
 * Heuristic search algorithms for sequence alignment and pattern matching
 */
class HeuristicSearch {
public:
    /**
     * A* search for optimal alignment path
     */
    class AStar {
    public:
        struct Node {
            int i, j;  // Position in alignment matrix
            int g;     // Cost from start
            int h;     // Heuristic cost to goal
            int f;     // Total cost (g + h)
            
            Node(int i, int j, int g, int h) : i(i), j(j), g(g), h(h), f(g + h) {}
            
            bool operator>(const Node& other) const {
                return f > other.f;
            }
        };
        
        /**
         * Find optimal alignment using A*
         * @param seq1 First sequence
         * @param seq2 Second sequence
         * @param match_score Match score
         * @param mismatch_score Mismatch score
         * @param gap_penalty Gap penalty
         * @return Alignment score
         */
        static int align(const std::string& seq1,
                        const std::string& seq2,
                        int match_score = 2,
                        int mismatch_score = -1,
                        int gap_penalty = -1);
        
    private:
        /**
         * Heuristic function (Manhattan distance approximation)
         */
        static int heuristic(int i, int j, int m, int n);
    };
    
    /**
     * Beam Search for sequence alignment
     */
    class BeamSearch {
    public:
        /**
         * Align sequences using beam search
         * @param seq1 First sequence
         * @param seq2 Second sequence
         * @param beam_width Beam width
         * @param match_score Match score
         * @param mismatch_score Mismatch score
         * @param gap_penalty Gap penalty
         * @return Alignment score
         */
        static int align(const std::string& seq1,
                        const std::string& seq2,
                        int beam_width = 10,
                        int match_score = 2,
                        int mismatch_score = -1,
                        int gap_penalty = -1);
    };
    
    /**
     * Genetic Algorithm for pattern evolution
     */
    class GeneticAlgorithm {
    public:
        struct Individual {
            std::string pattern;
            int fitness;
            
            Individual(const std::string& p, int f) : pattern(p), fitness(f) {}
        };
        
        /**
         * Evolve pattern to match sequence
         * @param sequence Target sequence
         * @param initial_pattern Initial pattern
         * @param population_size Population size
         * @param generations Number of generations
         * @param mutation_rate Mutation rate
         * @return Best evolved pattern
         */
        static std::string evolvePattern(const std::string& sequence,
                                        const std::string& initial_pattern,
                                        int population_size = 50,
                                        int generations = 100,
                                        double mutation_rate = 0.1);
        
    private:
        /**
         * Calculate fitness
         */
        static int calculateFitness(const std::string& sequence, const std::string& pattern);
        
        /**
         * Crossover two individuals
         */
        static Individual crossover(const Individual& parent1, const Individual& parent2);
        
        /**
         * Mutate individual
         */
        static Individual mutate(const Individual& individual, double mutation_rate);
        
        /**
         * Select parent (tournament selection)
         */
        static Individual tournamentSelection(const std::vector<Individual>& population);
    };
    
    /**
     * Simulated Annealing for alignment
     */
    class SimulatedAnnealing {
    public:
        /**
         * Align sequences using simulated annealing
         * @param seq1 First sequence
         * @param seq2 Second sequence
         * @param initial_temp Initial temperature
         * @param cooling_rate Cooling rate
         * @param iterations Number of iterations
         * @return Alignment score
         */
        static int align(const std::string& seq1,
                        const std::string& seq2,
                        double initial_temp = 100.0,
                        double cooling_rate = 0.95,
                        int iterations = 1000);
        
    private:
        /**
         * Generate neighbor solution
         */
        static std::pair<int, int> generateNeighbor(int i, int j, int m, int n);
        
        /**
         * Calculate alignment score at position
         */
        static int calculateScore(const std::string& seq1, const std::string& seq2,
                                 int i, int j, int match_score, int mismatch_score, int gap_penalty);
    };
};

#endif // HEURISTIC_SEARCH_H

