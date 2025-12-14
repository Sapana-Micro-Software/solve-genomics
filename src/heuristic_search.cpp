#include "heuristic_search.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <cctype>
#include <climits>

// A* Implementation

int HeuristicSearch::AStar::heuristic(int i, int j, int m, int n) {
    // Manhattan distance approximation
    return (m - i) + (n - j);
}

int HeuristicSearch::AStar::align(const std::string& seq1,
                                  const std::string& seq2,
                                  int match_score,
                                  int mismatch_score,
                                  int gap_penalty) {
    int m = static_cast<int>(seq1.length());
    int n = static_cast<int>(seq2.length());
    
    if (m == 0) return n * gap_penalty;
    if (n == 0) return m * gap_penalty;
    
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> open_set;
    std::map<std::pair<int, int>, int> g_score;
    
    Node start(0, 0, 0, heuristic(0, 0, m, n));
    open_set.push(start);
    g_score[{0, 0}] = 0;
    
    int best_score = INT_MIN;
    
    while (!open_set.empty()) {
        Node current = open_set.top();
        open_set.pop();
        
        if (current.i == m && current.j == n) {
            best_score = std::max(best_score, current.g);
            continue;
        }
        
        // Generate neighbors
        std::vector<Node> neighbors;
        
        // Match/Mismatch
        if (current.i < m && current.j < n) {
            int cost = (std::toupper(seq1[current.i]) == std::toupper(seq2[current.j])) 
                      ? match_score : mismatch_score;
            int new_g = current.g + cost;
            int h = heuristic(current.i + 1, current.j + 1, m, n);
            neighbors.push_back(Node(current.i + 1, current.j + 1, new_g, h));
        }
        
        // Gap in seq1
        if (current.i < m) {
            int new_g = current.g + gap_penalty;
            int h = heuristic(current.i + 1, current.j, m, n);
            neighbors.push_back(Node(current.i + 1, current.j, new_g, h));
        }
        
        // Gap in seq2
        if (current.j < n) {
            int new_g = current.g + gap_penalty;
            int h = heuristic(current.i, current.j + 1, m, n);
            neighbors.push_back(Node(current.i, current.j + 1, new_g, h));
        }
        
        for (const Node& neighbor : neighbors) {
            auto key = std::make_pair(neighbor.i, neighbor.j);
            if (g_score.find(key) == g_score.end() || g_score[key] < neighbor.g) {
                g_score[key] = neighbor.g;
                open_set.push(neighbor);
            }
        }
        
        // Limit search space
        if (open_set.size() > 10000) {
            break;
        }
    }
    
    return best_score;
}

// Beam Search Implementation

int HeuristicSearch::BeamSearch::align(const std::string& seq1,
                                       const std::string& seq2,
                                       int beam_width,
                                       int match_score,
                                       int mismatch_score,
                                       int gap_penalty) {
    int m = static_cast<int>(seq1.length());
    int n = static_cast<int>(seq2.length());
    
    if (m == 0) return n * gap_penalty;
    if (n == 0) return m * gap_penalty;
    
    // Beam: (i, j, score)
    std::vector<std::tuple<int, int, int>> beam;
    beam.push_back({0, 0, 0});
    
    int best_score = INT_MIN;
    
    for (int step = 0; step < m + n && !beam.empty(); ++step) {
        std::vector<std::tuple<int, int, int>> next_beam;
        std::map<std::pair<int, int>, int> best_at_pos;
        
        for (const auto& [i, j, score] : beam) {
            if (i == m && j == n) {
                best_score = std::max(best_score, score);
                continue;
            }
            
            // Match/Mismatch
            if (i < m && j < n) {
                int cost = (std::toupper(seq1[i]) == std::toupper(seq2[j])) 
                          ? match_score : mismatch_score;
                int new_score = score + cost;
                auto key = std::make_pair(i + 1, j + 1);
                if (best_at_pos.find(key) == best_at_pos.end() || 
                    best_at_pos[key] < new_score) {
                    best_at_pos[key] = new_score;
                    next_beam.push_back({i + 1, j + 1, new_score});
                }
            }
            
            // Gap in seq1
            if (i < m) {
                int new_score = score + gap_penalty;
                auto key = std::make_pair(i + 1, j);
                if (best_at_pos.find(key) == best_at_pos.end() || 
                    best_at_pos[key] < new_score) {
                    best_at_pos[key] = new_score;
                    next_beam.push_back({i + 1, j, new_score});
                }
            }
            
            // Gap in seq2
            if (j < n) {
                int new_score = score + gap_penalty;
                auto key = std::make_pair(i, j + 1);
                if (best_at_pos.find(key) == best_at_pos.end() || 
                    best_at_pos[key] < new_score) {
                    best_at_pos[key] = new_score;
                    next_beam.push_back({i, j + 1, new_score});
                }
            }
        }
        
        // Keep top beam_width
        std::sort(next_beam.begin(), next_beam.end(),
            [](const std::tuple<int, int, int>& a, const std::tuple<int, int, int>& b) {
                return std::get<2>(a) > std::get<2>(b);
            });
        
        if (static_cast<int>(next_beam.size()) > beam_width) {
            next_beam.resize(beam_width);
        }
        
        beam = next_beam;
    }
    
    return best_score;
}

// Genetic Algorithm Implementation

int HeuristicSearch::GeneticAlgorithm::calculateFitness(const std::string& sequence,
                                                         const std::string& pattern) {
    if (pattern.empty() || pattern.length() > sequence.length()) {
        return 0;
    }
    
    int matches = 0;
    for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {
        bool match = true;
        for (size_t j = 0; j < pattern.length(); ++j) {
            if (std::toupper(sequence[i + j]) != std::toupper(pattern[j])) {
                match = false;
                break;
            }
        }
        if (match) matches++;
    }
    
    return matches;
}

HeuristicSearch::GeneticAlgorithm::Individual 
HeuristicSearch::GeneticAlgorithm::crossover(const Individual& parent1, const Individual& parent2) {
    std::string child_pattern;
    size_t min_len = std::min(parent1.pattern.length(), parent2.pattern.length());
    
    for (size_t i = 0; i < min_len; ++i) {
        child_pattern += (i % 2 == 0) ? parent1.pattern[i] : parent2.pattern[i];
    }
    
    if (parent1.pattern.length() > min_len) {
        child_pattern += parent1.pattern.substr(min_len);
    } else if (parent2.pattern.length() > min_len) {
        child_pattern += parent2.pattern.substr(min_len);
    }
    
    return Individual(child_pattern, 0);
}

HeuristicSearch::GeneticAlgorithm::Individual 
HeuristicSearch::GeneticAlgorithm::mutate(const Individual& individual, double mutation_rate) {
    std::string mutated = individual.pattern;
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::uniform_int_distribution<int> base_dist(0, 3);
    
    const char bases[] = {'A', 'T', 'G', 'C'};
    
    for (char& c : mutated) {
        if (dist(rng) < mutation_rate) {
            c = bases[base_dist(rng)];
        }
    }
    
    return Individual(mutated, 0);
}

HeuristicSearch::GeneticAlgorithm::Individual 
HeuristicSearch::GeneticAlgorithm::tournamentSelection(const std::vector<Individual>& population) {
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, population.size() - 1);
    
    const Individual* best = &population[dist(rng)];
    for (int i = 0; i < 3; ++i) {
        const Individual* candidate = &population[dist(rng)];
        if (candidate->fitness > best->fitness) {
            best = candidate;
        }
    }
    
    return *best;
}

std::string HeuristicSearch::GeneticAlgorithm::evolvePattern(const std::string& sequence,
                                                            const std::string& initial_pattern,
                                                            int population_size,
                                                            int generations,
                                                            double mutation_rate) {
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> len_dist(initial_pattern.length() - 2, 
                                                   initial_pattern.length() + 2);
    
    // Initialize population
    std::vector<Individual> population;
    population.push_back(Individual(initial_pattern, calculateFitness(sequence, initial_pattern)));
    
    for (int i = 1; i < population_size; ++i) {
        size_t len = len_dist(rng);
        std::string pattern;
        const char bases[] = {'A', 'T', 'G', 'C'};
        std::uniform_int_distribution<int> base_dist(0, 3);
        
        for (size_t j = 0; j < len; ++j) {
            pattern += bases[base_dist(rng)];
        }
        
        int fitness = calculateFitness(sequence, pattern);
        population.push_back(Individual(pattern, fitness));
    }
    
    // Evolve
    for (int gen = 0; gen < generations; ++gen) {
        // Evaluate fitness
        for (auto& individual : population) {
            individual.fitness = calculateFitness(sequence, individual.pattern);
        }
        
        // Sort by fitness
        std::sort(population.begin(), population.end(),
            [](const Individual& a, const Individual& b) {
                return a.fitness > b.fitness;
            });
        
        // Create new generation
        std::vector<Individual> new_population;
        
        // Keep top 10%
        int elite_size = population_size / 10;
        for (int i = 0; i < elite_size; ++i) {
            new_population.push_back(population[i]);
        }
        
        // Generate offspring
        while (static_cast<int>(new_population.size()) < population_size) {
            Individual parent1 = tournamentSelection(population);
            Individual parent2 = tournamentSelection(population);
            Individual child = crossover(parent1, parent2);
            child = mutate(child, mutation_rate);
            child.fitness = calculateFitness(sequence, child.pattern);
            new_population.push_back(child);
        }
        
        population = new_population;
    }
    
    // Return best
    std::sort(population.begin(), population.end(),
        [](const Individual& a, const Individual& b) {
            return a.fitness > b.fitness;
        });
    
    return population[0].pattern;
}

// Simulated Annealing Implementation

std::pair<int, int> HeuristicSearch::SimulatedAnnealing::generateNeighbor(int i, int j, int m, int n) {
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> move_dist(0, 2);
    
    int move = move_dist(rng);
    int new_i = i, new_j = j;
    
    switch (move) {
        case 0:  // Move diagonal
            if (i < m && j < n) {
                new_i++; new_j++;
            }
            break;
        case 1:  // Move down
            if (i < m) new_i++;
            break;
        case 2:  // Move right
            if (j < n) new_j++;
            break;
    }
    
    return {new_i, new_j};
}

int HeuristicSearch::SimulatedAnnealing::calculateScore(const std::string& seq1,
                                                        const std::string& seq2,
                                                        int i, int j,
                                                        int match_score,
                                                        int mismatch_score,
                                                        int gap_penalty) {
    int score = 0;
    
    // Simplified: calculate score up to position (i, j)
    int pos1 = 0, pos2 = 0;
    
    while (pos1 < i && pos2 < j) {
        if (std::toupper(seq1[pos1]) == std::toupper(seq2[pos2])) {
            score += match_score;
        } else {
            score += mismatch_score;
        }
        pos1++; pos2++;
    }
    
    score += (i - pos1) * gap_penalty;
    score += (j - pos2) * gap_penalty;
    
    return score;
}

int HeuristicSearch::SimulatedAnnealing::align(const std::string& seq1,
                                               const std::string& seq2,
                                               double initial_temp,
                                               double cooling_rate,
                                               int iterations) {
    int m = static_cast<int>(seq1.length());
    int n = static_cast<int>(seq2.length());
    
    int match_score = 2, mismatch_score = -1, gap_penalty = -1;
    
    int current_i = 0, current_j = 0;
    int current_score = 0;
    int best_score = current_score;
    
    double temperature = initial_temp;
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    
    for (int iter = 0; iter < iterations; ++iter) {
        auto [new_i, new_j] = generateNeighbor(current_i, current_j, m, n);
        
        int new_score = calculateScore(seq1, seq2, new_i, new_j, 
                                      match_score, mismatch_score, gap_penalty);
        
        int delta = new_score - current_score;
        
        if (delta > 0 || prob_dist(rng) < std::exp(delta / temperature)) {
            current_i = new_i;
            current_j = new_j;
            current_score = new_score;
            
            if (current_score > best_score) {
                best_score = current_score;
            }
        }
        
        temperature *= cooling_rate;
        
        if (current_i == m && current_j == n) {
            break;
        }
    }
    
    return best_score;
}

