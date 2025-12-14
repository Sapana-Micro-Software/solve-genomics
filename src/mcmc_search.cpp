#include "mcmc_search.h"
#include "ExactMatch.h"
#include <algorithm>
#include <cmath>
#include <random>

MCMCSearch::MCMCSearch(int max_iterations, double mutation_rate, double temperature)
    : max_iterations_(max_iterations), 
      mutation_rate_(mutation_rate),
      temperature_(temperature),
      rng_(std::random_device{}()),
      uniform_dist_(0.0, 1.0) {
}

MCMCSearch::MCMCResult MCMCSearch::search(const std::string& sequence, const std::string& initial_pattern) {
    MCMCResult result;
    
    if (sequence.empty() || initial_pattern.empty()) {
        return result;
    }
    
    std::string current_pattern = initial_pattern;
    int current_fitness = calculateFitness(sequence, current_pattern);
    
    std::string best_pattern = current_pattern;
    int best_fitness = current_fitness;
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        // Propose new pattern by mutation
        std::string proposed_pattern = mutatePattern(current_pattern, mutation_rate_);
        int proposed_fitness = calculateFitness(sequence, proposed_pattern);
        
        // Metropolis-Hastings acceptance
        if (acceptProposal(current_fitness, proposed_fitness, temperature_)) {
            current_pattern = proposed_pattern;
            current_fitness = proposed_fitness;
            
            // Update best if better
            if (current_fitness > best_fitness) {
                best_pattern = current_pattern;
                best_fitness = current_fitness;
            }
        }
        
        // If we found exact matches, we can stop early
        if (best_fitness > 0 && current_pattern == initial_pattern) {
            // Check if evolved pattern matches exactly
            ExactMatch matcher;
            SearchResult exact_result = matcher.search(sequence, best_pattern);
            if (exact_result.count > 0) {
                result.positions.clear();
                for (int pos : exact_result.positions) {
                    result.positions.push_back(static_cast<size_t>(pos));
                }
                result.count = exact_result.count;
                result.evolved_pattern = best_pattern;
                result.iterations = iter + 1;
                result.final_fitness = best_fitness;
                break;
            }
        }
    }
    
    // Final search with best pattern
    ExactMatch matcher;
    SearchResult final_result = matcher.search(sequence, best_pattern);
    result.positions.clear();
    for (int pos : final_result.positions) {
        result.positions.push_back(static_cast<size_t>(pos));
    }
    result.count = final_result.count;
    result.evolved_pattern = best_pattern;
    result.iterations = max_iterations_;
    result.final_fitness = best_fitness;
    
    return result;
}

int MCMCSearch::calculateFitness(const std::string& sequence, const std::string& pattern) {
    if (pattern.empty() || pattern.length() > sequence.length()) {
        return 0;
    }
    
    ExactMatch matcher;
    SearchResult result = matcher.search(sequence, pattern);
    return result.count;
}

std::string MCMCSearch::mutatePattern(const std::string& pattern, double mutation_rate) {
    std::string mutated = pattern;
    
    for (size_t i = 0; i < mutated.length(); ++i) {
        if (uniform_dist_(rng_) < mutation_rate) {
            // Apply random mutation
            MutationType type = static_cast<MutationType>(
                std::uniform_int_distribution<int>(0, 2)(rng_));
            mutated = applyMutation(mutated, type, i);
        }
    }
    
    return mutated;
}

bool MCMCSearch::acceptProposal(int current_fitness, int proposed_fitness, double temperature) {
    if (proposed_fitness > current_fitness) {
        return true;  // Always accept better proposals
    }
    
    // Accept worse proposals with probability based on temperature
    double delta = proposed_fitness - current_fitness;
    double probability = std::exp(delta / temperature);
    
    return uniform_dist_(rng_) < probability;
}

MCMCSearch::MCMCResult MCMCSearch::simulatedAnnealingSearch(const std::string& sequence,
                                                            const std::string& initial_pattern) {
    MCMCResult result;
    
    if (sequence.empty() || initial_pattern.empty()) {
        return result;
    }
    
    std::string current_pattern = initial_pattern;
    int current_fitness = calculateFitness(sequence, current_pattern);
    
    std::string best_pattern = current_pattern;
    int best_fitness = current_fitness;
    
    double initial_temp = temperature_;
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        // Cool temperature
        double current_temp = coolTemperature(iter, max_iterations_);
        
        // Propose mutation
        std::string proposed_pattern = mutatePattern(current_pattern, mutation_rate_);
        int proposed_fitness = calculateFitness(sequence, proposed_pattern);
        
        // Accept with current temperature
        if (acceptProposal(current_fitness, proposed_fitness, current_temp)) {
            current_pattern = proposed_pattern;
            current_fitness = proposed_fitness;
            
            if (current_fitness > best_fitness) {
                best_pattern = current_pattern;
                best_fitness = current_fitness;
            }
        }
    }
    
    // Final search
    ExactMatch matcher;
    SearchResult final_result = matcher.search(sequence, best_pattern);
    result.positions.clear();
    for (int pos : final_result.positions) {
        result.positions.push_back(static_cast<size_t>(pos));
    }
    result.count = final_result.count;
    result.evolved_pattern = best_pattern;
    result.iterations = max_iterations_;
    result.final_fitness = best_fitness;
    
    return result;
}

std::string MCMCSearch::applyMutation(const std::string& pattern, MutationType type, size_t position) {
    std::string mutated = pattern;
    
    if (position >= mutated.length()) {
        return mutated;
    }
    
    switch (type) {
        case SUBSTITUTION: {
            char current = std::toupper(mutated[position]);
            char new_base = randomBase();
            // Ensure it's different
            while (new_base == current) {
                new_base = randomBase();
            }
            mutated[position] = new_base;
            break;
        }
        case INSERTION: {
            char new_base = randomBase();
            mutated.insert(position, 1, new_base);
            break;
        }
        case DELETION: {
            if (mutated.length() > 1) {
                mutated.erase(position, 1);
            }
            break;
        }
    }
    
    return mutated;
}

char MCMCSearch::randomBase() {
    const char bases[] = {'A', 'T', 'G', 'C'};
    std::uniform_int_distribution<int> base_dist(0, 3);
    return bases[base_dist(rng_)];
}

double MCMCSearch::coolTemperature(int iteration, int max_iter) {
    // Exponential cooling schedule
    return temperature_ * std::exp(-5.0 * iteration / max_iter);
}

