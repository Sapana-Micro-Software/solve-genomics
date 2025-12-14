#ifndef PROBABILISTIC_ML_H
#define PROBABILISTIC_ML_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>
#include <random>
#include <cmath>

/**
 * Probabilistic Machine Learning Methods for DNA Sequence Pattern Matching
 * Includes: Naive Bayes, Deep Belief Networks (DBN), Restricted Boltzmann Machines (RBM),
 *           Markov Decision Processes (MDP), and Markov Random Fields (MRF)
 */
namespace ProbabilisticML {

    /**
     * Naive Bayes Classifier for Sequence Pattern Matching
     * Assumes independence between sequence positions (naive assumption)
     */
    class NaiveBayes {
    public:
        struct ClassificationResult {
            std::string predicted_class;                                                      // Predicted pattern class
            double confidence;                                                                // Confidence score (0-1)
            std::map<std::string, double> class_probabilities;                                // Probability for each class
            
            ClassificationResult() : confidence(0.0) {}
        };
        
        NaiveBayes();
        
        /**
         * Train classifier on labeled sequences
         * @param sequences Training sequences
         * @param labels Class labels for each sequence
         */
        void train(const std::vector<std::string>& sequences,
                  const std::vector<std::string>& labels);
        
        /**
         * Classify a sequence (find matching pattern class)
         * @param sequence Sequence to classify
         * @return ClassificationResult with predicted class and confidence
         */
        ClassificationResult classify(const std::string& sequence);
        
        /**
         * Calculate probability of sequence given class
         * @param sequence Input sequence
         * @param class_label Class label
         * @return Probability P(sequence | class)
         */
        double probabilityGivenClass(const std::string& sequence, const std::string& class_label);
        
    private:
        std::map<std::string, int> class_counts_;                                             // Count of each class
        std::map<std::string, std::map<char, int>> class_base_counts_;                        // Base counts per class
        std::map<std::string, std::map<int, std::map<char, double>>> position_probs_;         // Position-specific probabilities
        int total_sequences_;                                                                 // Total training sequences
        double smoothing_alpha_;                                                               // Laplace smoothing parameter
        
        /**
         * Calculate base probability with smoothing
         */
        double getBaseProbability(const std::string& class_label, int position, char base);
    };
    
    /**
     * Restricted Boltzmann Machine (RBM) for Sequence Feature Learning
     * Unsupervised learning of sequence representations
     */
    class RestrictedBoltzmannMachine {
    public:
        struct RBMParameters {
            int visible_size;                                                                 // Size of visible layer (sequence encoding)
            int hidden_size;                                                                  // Size of hidden layer (features)
            double learning_rate;                                                             // Learning rate for training
            int num_epochs;                                                                   // Number of training epochs
            
            RBMParameters() : visible_size(64), hidden_size(32), learning_rate(0.1), num_epochs(100) {}
        };
        
        RestrictedBoltzmannMachine(int visible_size, int hidden_size);
        
        /**
         * Train RBM on sequences (unsupervised)
         * @param sequences Training sequences
         */
        void train(const std::vector<std::string>& sequences);
        
        /**
         * Encode sequence to hidden representation
         * @param sequence Input sequence
         * @return Hidden layer activations (feature vector)
         */
        std::vector<double> encode(const std::string& sequence);
        
        /**
         * Reconstruct sequence from hidden representation
         * @param hidden Hidden layer activations
         * @return Reconstructed sequence
         */
        std::string decode(const std::vector<double>& hidden);
        
        /**
         * Find similar sequences using RBM features
         * @param query_sequence Query sequence
         * @param database_sequences Database of sequences
         * @param top_k Number of similar sequences to return
         * @return Vector of (index, similarity_score) pairs
         */
        std::vector<std::pair<size_t, double>> findSimilar(const std::string& query_sequence,
                                                           const std::vector<std::string>& database_sequences,
                                                           int top_k = 5);
        
    private:
        int visible_size_;                                                                    // Visible layer size
        int hidden_size_;                                                                     // Hidden layer size
        std::vector<std::vector<double>> weights_;                                            // Weight matrix (visible x hidden)
        std::vector<double> visible_bias_;                                                    // Visible layer biases
        std::vector<double> hidden_bias_;                                                     // Hidden layer biases
        double learning_rate_;                                                                // Learning rate
        std::mt19937 rng_;                                                                    // Random number generator
        
        /**
         * Encode sequence to vector
         */
        std::vector<double> sequenceToVector(const std::string& sequence);
        
        /**
         * Vector to sequence
         */
        std::string vectorToSequence(const std::vector<double>& vector);
        
        /**
         * Sigmoid activation function
         */
        double sigmoid(double x);
        
        /**
         * Sample from Bernoulli distribution
         */
        double sampleBernoulli(double p);
        
        /**
         * Contrastive Divergence (CD-k) training step
         */
        void contrastiveDivergence(const std::vector<double>& visible, int k = 1);
    };
    
    /**
     * Deep Belief Network (DBN) for Hierarchical Feature Learning
     * Stacked RBMs for deep representation learning
     */
    class DeepBeliefNetwork {
    public:
        struct DBNParameters {
            std::vector<int> layer_sizes;                                                      // Size of each hidden layer
            double learning_rate;                                                             // Learning rate
            int num_epochs_per_layer;                                                         // Epochs for each RBM layer
            
            DBNParameters() : learning_rate(0.1), num_epochs_per_layer(50) {}
        };
        
        DeepBeliefNetwork(const std::vector<int>& layer_sizes);
        
        /**
         * Train DBN layer by layer (greedy training)
         * @param sequences Training sequences
         */
        void train(const std::vector<std::string>& sequences);
        
        /**
         * Encode sequence through all layers
         * @param sequence Input sequence
         * @return Final hidden representation
         */
        std::vector<double> encode(const std::string& sequence);
        
        /**
         * Find patterns using deep features
         * @param sequence Query sequence
         * @param pattern_patterns Pattern database
         * @return Similarity scores for each pattern
         */
        std::vector<std::pair<std::string, double>> matchPatterns(const std::string& sequence,
                                                                  const std::vector<std::string>& pattern_patterns);
        
    private:
        std::vector<int> layer_sizes_;                                                        // Size of each layer
        std::vector<RestrictedBoltzmannMachine> rbms_;                                         // Stack of RBMs
        int visible_size_;                                                                    // Input size
        
        /**
         * Encode sequence to vector
         */
        std::vector<double> sequenceToVector(const std::string& sequence);
    };
    
    /**
     * Markov Decision Process (MDP) for Sequential Pattern Matching
     * Models pattern matching as a decision process with states and actions
     */
    class MarkovDecisionProcess {
    public:
        struct State {
            int position;                                                                     // Current position in sequence
            int pattern_index;                                                                // Current position in pattern
            int matches;                                                                      // Number of matches so far
            int errors;                                                                       // Number of errors allowed
            
            bool operator<(const State& other) const {
                if (position != other.position) return position < other.position;
                if (pattern_index != other.pattern_index) return pattern_index < other.pattern_index;
                if (matches != other.matches) return matches < other.matches;
                return errors < other.errors;
            }
        };
        
        struct Action {
            enum Type { MATCH, INSERT, DELETE, SKIP };
            Type type;                                                                        // Action type
            double reward;                                                                    // Reward for this action
        };
        
        struct MDPResult {
            std::vector<int> positions;                                                       // Matched positions
            double total_reward;                                                              // Total reward
            std::vector<Action> action_sequence;                                             // Sequence of actions taken
            
            MDPResult() : total_reward(0.0) {}
        };
        
        MarkovDecisionProcess(double match_reward = 2.0, double mismatch_penalty = -1.0,
                             double gap_penalty = -1.0);
        
        /**
         * Find pattern using MDP policy
         * @param sequence Sequence to search in
         * @param pattern Pattern to find
         * @param max_errors Maximum errors allowed
         * @return MDPResult with matched positions and reward
         */
        MDPResult findPattern(const std::string& sequence, const std::string& pattern, int max_errors = 0);
        
        /**
         * Value iteration for optimal policy
         * @param sequence Sequence
         * @param pattern Pattern
         * @return Optimal value function
         */
        std::map<State, double> valueIteration(const std::string& sequence, const std::string& pattern);
        
    private:
        double match_reward_;                                                                 // Reward for match
        double mismatch_penalty_;                                                             // Penalty for mismatch
        double gap_penalty_;                                                                  // Penalty for gap
        double discount_factor_;                                                              // Discount factor for future rewards
        
        /**
         * Get reward for action in state
         */
        double getReward(const State& state, const Action& action, char seq_char, char pat_char);
        
        /**
         * Get next state from current state and action
         */
        State getNextState(const State& current, const Action& action);
        
        /**
         * Get possible actions from state
         */
        std::vector<Action> getPossibleActions(const State& state, const std::string& sequence, const std::string& pattern);
    };
    
    /**
     * Markov Random Field (MRF) for Sequence Pattern Matching
     * Models spatial dependencies between sequence positions
     */
    class MarkovRandomField {
    public:
        struct MRFResult {
            std::vector<int> matched_positions;                                               // Positions where pattern matches
            double energy;                                                                    // Energy of configuration
            std::map<int, double> position_probabilities;                                     // Probability of match at each position
            
            MRFResult() : energy(0.0) {}
        };
        
        MarkovRandomField(double match_energy = -2.0, double mismatch_energy = 1.0,
                         double smoothness_weight = 0.5);
        
        /**
         * Find pattern using MRF inference
         * @param sequence Sequence to search in
         * @param pattern Pattern to find
         * @return MRFResult with matched positions and probabilities
         */
        MRFResult findPattern(const std::string& sequence, const std::string& pattern);
        
        /**
         * Gibbs sampling for MRF inference
         * @param sequence Sequence
         * @param pattern Pattern
         * @param num_iterations Number of sampling iterations
         * @return MRFResult with sampled configuration
         */
        MRFResult gibbsSampling(const std::string& sequence, const std::string& pattern, int num_iterations = 1000);
        
        /**
         * Belief propagation for MRF inference
         * @param sequence Sequence
         * @param pattern Pattern
         * @return MRFResult with inferred probabilities
         */
        MRFResult beliefPropagation(const std::string& sequence, const std::string& pattern);
        
    private:
        double match_energy_;                                                                 // Energy for match
        double mismatch_energy_;                                                              // Energy for mismatch
        double smoothness_weight_;                                                            // Weight for spatial smoothness
        std::mt19937 rng_;                                                                    // Random number generator
        
        /**
         * Calculate energy of configuration
         */
        double calculateEnergy(const std::vector<bool>& matches, const std::string& sequence, const std::string& pattern);
        
        /**
         * Calculate local energy at position
         */
        double localEnergy(int position, bool is_match, const std::vector<bool>& neighbors,
                          const std::string& sequence, const std::string& pattern);
        
        /**
         * Sample from conditional distribution
         */
        bool sampleConditional(int position, const std::vector<bool>& current_state,
                              const std::string& sequence, const std::string& pattern);
    };
    
} // namespace ProbabilisticML

#endif // PROBABILISTIC_ML_H

