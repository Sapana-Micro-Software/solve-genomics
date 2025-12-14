#include "probabilistic_ml.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cctype>
#include <limits>

namespace ProbabilisticML {

    // ==================== Naive Bayes Implementation ====================
    
    NaiveBayes::NaiveBayes() : total_sequences_(0), smoothing_alpha_(1.0) {
    }
    
    void NaiveBayes::train(const std::vector<std::string>& sequences,
                          const std::vector<std::string>& labels) {
        if (sequences.size() != labels.size() || sequences.empty()) {                        // Validate input: check sizes match and not empty
            return;                                                                            // Return early if invalid input
        }
        
        total_sequences_ = static_cast<int>(sequences.size());                                 // Store total number of training sequences
        class_counts_.clear();                                                                 // Clear previous class frequency counts
        class_base_counts_.clear();                                                            // Clear previous base frequency counts per class
        position_probs_.clear();                                                               // Clear previous position-specific probabilities
        
        // Count classes and bases
        for (size_t i = 0; i < sequences.size(); ++i) {                                       // Iterate through all training sequences
            const std::string& seq = sequences[i];                                            // Get current training sequence
            const std::string& label = labels[i];                                            // Get corresponding class label
            
            class_counts_[label]++;                                                           // Increment count for this class (prior probability)
            
            if (class_base_counts_.find(label) == class_base_counts_.end()) {                 // Check if this is a new class
                class_base_counts_[label] = std::map<char, int>();                            // Initialize base count map for new class
            }
            
            // Count bases at each position
            for (size_t pos = 0; pos < seq.length(); ++pos) {                                 // Iterate through each position in sequence
                char base = std::toupper(seq[pos]);                                           // Get base at current position (uppercase)
                class_base_counts_[label][base]++;                                            // Increment count of this base for this class
                
                if (position_probs_[label].find(static_cast<int>(pos)) == position_probs_[label].end()) { // Check if position initialized
                    position_probs_[label][static_cast<int>(pos)] = std::map<char, double>(); // Initialize probability map for this position
                }
            }
        }
        
        // Calculate probabilities with Laplace smoothing
        for (const auto& class_pair : class_counts_) {                                        // Iterate through each class in training data
            const std::string& label = class_pair.first;                                      // Get class label
            int class_count = class_pair.second;                                               // Get number of sequences in this class
            
            for (auto& pos_pair : position_probs_[label]) {                                    // Iterate through each position in sequences
                int pos = pos_pair.first;                                                     // Get current position index
                double total_bases = 0.0;                                                      // Initialize total base count (for normalization)
                
                // Count bases at this position for this class
                for (char base : {'A', 'T', 'G', 'C'}) {                                      // Iterate through all four DNA bases
                    int count = 0;                                                             // Initialize count for this base at this position
                    for (size_t i = 0; i < sequences.size(); ++i) {                           // Count occurrences in training sequences
                        if (labels[i] == label && pos < static_cast<int>(sequences[i].length())) { // Check if sequence belongs to class and has this position
                            if (std::toupper(sequences[i][pos]) == base) {                     // Check if base at position matches
                                count++;                                                       // Increment count for this base
                            }
                        }
                    }
                    total_bases += count;                                                      // Add to total base count
                    pos_pair.second[base] = (count + smoothing_alpha_) / (class_count + 4 * smoothing_alpha_); // Calculate probability with Laplace smoothing
                }
            }
        }
    }
    
    NaiveBayes::ClassificationResult NaiveBayes::classify(const std::string& sequence) {
        ClassificationResult result;                                                          // Initialize classification result structure
        
        if (class_counts_.empty() || sequence.empty()) {                                      // Check if classifier is trained and input is valid
            return result;                                                                    // Return empty result if not trained or invalid input
        }
        
        double total_prob = 0.0;                                                              // Initialize total probability (for normalization)
        std::map<std::string, double> class_log_probs;                                        // Store log probabilities for each class (for numerical stability)
        
        // Calculate log probability for each class
        for (const auto& class_pair : class_counts_) {                                        // Iterate through each class in training data
            const std::string& label = class_pair.first;                                      // Get class label
            double class_prior = static_cast<double>(class_pair.second) / total_sequences_;   // Calculate prior probability P(class) = count(class) / total
            double log_prob = std::log(class_prior);                                          // Start with log of prior probability (Bayes' theorem)
            
            // Calculate likelihood P(sequence | class) using naive assumption (independence)
            for (size_t pos = 0; pos < sequence.length(); ++pos) {                            // Iterate through each position in query sequence
                char base = std::toupper(sequence[pos]);                                      // Get base at current position (uppercase)
                double base_prob = getBaseProbability(label, static_cast<int>(pos), base);    // Get probability P(base | class, position)
                if (base_prob > 0.0) {                                                        // Check if probability is non-zero (avoid log(0))
                    log_prob += std::log(base_prob);                                          // Add log probability (product becomes sum in log space)
                }
            }
            
            class_log_probs[label] = log_prob;                                                // Store log probability for this class
        }
        
        // Find class with maximum probability
        double max_log_prob = -std::numeric_limits<double>::infinity();                       // Initialize maximum log probability to negative infinity
        for (const auto& prob_pair : class_log_probs) {                                       // Iterate through all class log probabilities
            if (prob_pair.second > max_log_prob) {                                            // Check if this class has higher log probability
                max_log_prob = prob_pair.second;                                              // Update maximum log probability
                result.predicted_class = prob_pair.first;                                     // Update predicted class to this class
            }
        }
        
        // Convert log probabilities to probabilities and normalize
        double sum_exp = 0.0;                                                                 // Initialize sum of exponentials (for normalization)
        for (auto& prob_pair : class_log_probs) {                                             // Iterate through all class log probabilities
            double exp_prob = std::exp(prob_pair.second - max_log_prob);                      // Exponentiate with numerical stability (subtract max)
            prob_pair.second = exp_prob;                                                      // Store probability (unnormalized)
            sum_exp += exp_prob;                                                              // Add to sum for normalization
        }
        
        // Normalize probabilities
        for (auto& prob_pair : class_log_probs) {                                            // Iterate through all class probabilities
            prob_pair.second /= sum_exp;                                                      // Normalize by dividing by sum (Bayes' rule)
            result.class_probabilities[prob_pair.first] = prob_pair.second;                   // Store normalized probability in result
        }
        
        result.confidence = result.class_probabilities[result.predicted_class];                // Set confidence to probability of predicted class
        
        return result;                                                                        // Return complete classification result
    }
    
    double NaiveBayes::getBaseProbability(const std::string& class_label, int position, char base) {
        if (position_probs_.find(class_label) == position_probs_.end()) {                    // Check if class exists
            return 0.25;                                                                      // Return uniform prior if not found
        }
        
        if (position_probs_[class_label].find(position) == position_probs_[class_label].end()) { // Check if position exists
            return 0.25;                                                                      // Return uniform prior if not found
        }
        
        if (position_probs_[class_label][position].find(base) == position_probs_[class_label][position].end()) { // Check if base exists
            return 0.25;                                                                      // Return uniform prior if not found
        }
        
        return position_probs_[class_label][position][base];                                  // Return stored probability
    }
    
    // ==================== Restricted Boltzmann Machine Implementation ====================
    
    RestrictedBoltzmannMachine::RestrictedBoltzmannMachine(int visible_size, int hidden_size)
        : visible_size_(visible_size), hidden_size_(hidden_size), learning_rate_(0.1),
          rng_(std::random_device{}()) {
        // Initialize weights with small random values
        weights_.resize(visible_size_);                                                        // Resize weight matrix
        for (int i = 0; i < visible_size_; ++i) {                                            // Initialize each row
            weights_[i].resize(hidden_size_);                                                 // Resize row
            for (int j = 0; j < hidden_size_; ++j) {                                         // Initialize each weight
                std::uniform_real_distribution<double> dist(-0.1, 0.1);                       // Small random initialization
                weights_[i][j] = dist(rng_);                                                  // Set random weight
            }
        }
        
        visible_bias_.resize(visible_size_, 0.0);                                            // Initialize visible biases
        hidden_bias_.resize(hidden_size_, 0.0);                                               // Initialize hidden biases
    }
    
    void RestrictedBoltzmannMachine::train(const std::vector<std::string>& sequences) {
        if (sequences.empty()) {                                                               // Check if training sequences are provided
            return;                                                                            // Return early if no sequences
        }
        
        const int num_epochs = 100;                                                            // Number of training epochs (iterations over dataset)
        const int k = 1;                                                                      // CD-k parameter (k steps of Gibbs sampling)
        
        for (int epoch = 0; epoch < num_epochs; ++epoch) {                                    // Train for specified number of epochs
            for (const auto& sequence : sequences) {                                           // Iterate through each training sequence
                std::vector<double> visible = sequenceToVector(sequence);                     // Convert DNA sequence to numerical vector
                if (visible.size() != static_cast<size_t>(visible_size_)) {                   // Check if vector size matches expected visible size
                    continue;                                                                 // Skip this sequence if size mismatch
                }
                contrastiveDivergence(visible, k);                                            // Perform Contrastive Divergence update (CD-k)
            }
        }
    }
    
    std::vector<double> RestrictedBoltzmannMachine::encode(const std::string& sequence) {
        std::vector<double> visible = sequenceToVector(sequence);                             // Convert DNA sequence to numerical vector representation
        if (visible.size() != static_cast<size_t>(visible_size_)) {                           // Check if vector size matches expected visible layer size
            visible.resize(visible_size_, 0.0);                                               // Resize to expected size (pad with zeros if needed)
        }
        
        std::vector<double> hidden(hidden_size_, 0.0);                                        // Initialize hidden layer activations to zero
        
        // Calculate hidden activations
        for (int j = 0; j < hidden_size_; ++j) {                                             // Iterate through each hidden unit
            double activation = hidden_bias_[j];                                               // Start with hidden bias for this unit
            for (int i = 0; i < visible_size_; ++i) {                                        // Sum contributions from all visible units
                activation += visible[i] * weights_[i][j];                                    // Add weighted input: visible[i] * weight[i][j]
            }
            hidden[j] = sigmoid(activation);                                                  // Apply sigmoid activation function: 1/(1+exp(-x))
        }
        
        return hidden;                                                                        // Return hidden layer representation (feature vector)
    }
    
    std::string RestrictedBoltzmannMachine::decode(const std::vector<double>& hidden) {
        if (hidden.size() != static_cast<size_t>(hidden_size_)) {                             // Check size match
            return "";                                                                        // Return empty if mismatch
        }
        
        std::vector<double> visible(visible_size_, 0.0);                                      // Initialize visible layer
        
        // Calculate visible activations
        for (int i = 0; i < visible_size_; ++i) {                                            // For each visible unit
            double activation = visible_bias_[i];                                              // Start with bias
            for (int j = 0; j < hidden_size_; ++j) {                                         // Sum over hidden units
                activation += hidden[j] * weights_[i][j];                                     // Weighted sum
            }
            visible[i] = sigmoid(activation);                                                 // Apply sigmoid activation
        }
        
        return vectorToSequence(visible);                                                     // Convert back to sequence
    }
    
    std::vector<std::pair<size_t, double>> RestrictedBoltzmannMachine::findSimilar(
            const std::string& query_sequence,
            const std::vector<std::string>& database_sequences,
            int top_k) {
        std::vector<double> query_hidden = encode(query_sequence);                            // Encode query sequence
        
        std::vector<std::pair<size_t, double>> similarities;                                  // Store similarities
        for (size_t i = 0; i < database_sequences.size(); ++i) {                              // Iterate through database
            std::vector<double> db_hidden = encode(database_sequences[i]);                   // Encode database sequence
            
            // Calculate cosine similarity
            double dot_product = 0.0;                                                         // Initialize dot product
            double norm_query = 0.0;                                                          // Initialize query norm
            double norm_db = 0.0;                                                             // Initialize database norm
            
            for (size_t j = 0; j < query_hidden.size() && j < db_hidden.size(); ++j) {      // Calculate similarity
                dot_product += query_hidden[j] * db_hidden[j];                               // Dot product
                norm_query += query_hidden[j] * query_hidden[j];                              // Query norm squared
                norm_db += db_hidden[j] * db_hidden[j];                                      // Database norm squared
            }
            
            double similarity = 0.0;                                                          // Initialize similarity
            if (norm_query > 0.0 && norm_db > 0.0) {                                         // Avoid division by zero
                similarity = dot_product / (std::sqrt(norm_query) * std::sqrt(norm_db));     // Cosine similarity
            }
            
            similarities.push_back({i, similarity});                                          // Store similarity
        }
        
        // Sort by similarity (descending)
        std::sort(similarities.begin(), similarities.end(),
                  [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                      return a.second > b.second;                                              // Sort descending
                  });
        
        // Return top-k
        if (top_k > 0 && top_k < static_cast<int>(similarities.size())) {                    // Check if need to limit
            similarities.resize(top_k);                                                       // Resize to top-k
        }
        
        return similarities;                                                                  // Return similar sequences
    }
    
    std::vector<double> RestrictedBoltzmannMachine::sequenceToVector(const std::string& sequence) {
        std::vector<double> vector;                                                           // Initialize vector
        vector.reserve(sequence.length() * 4);                                               // Reserve space (4 bases)
        
        for (char c : sequence) {                                                             // Iterate through sequence
            char base = std::toupper(c);                                                      // Convert to uppercase
            // One-hot encoding
            vector.push_back(base == 'A' ? 1.0 : 0.0);                                       // A
            vector.push_back(base == 'T' ? 1.0 : 0.0);                                       // T
            vector.push_back(base == 'G' ? 1.0 : 0.0);                                       // G
            vector.push_back(base == 'C' ? 1.0 : 0.0);                                       // C
        }
        
        // Pad or truncate to visible_size_
        if (vector.size() < static_cast<size_t>(visible_size_)) {                            // Check if need padding
            vector.resize(visible_size_, 0.0);                                               // Pad with zeros
        } else if (vector.size() > static_cast<size_t>(visible_size_)) {                     // Check if need truncation
            vector.resize(visible_size_);                                                     // Truncate
        }
        
        return vector;                                                                        // Return encoded vector
    }
    
    std::string RestrictedBoltzmannMachine::vectorToSequence(const std::vector<double>& vector) {
        std::string sequence;                                                                 // Initialize sequence
        
        for (size_t i = 0; i < vector.size(); i += 4) {                                      // Process in groups of 4
            if (i + 3 >= vector.size()) break;                                                // Check bounds
            
            // Find maximum (one-hot decoding)
            int max_idx = 0;                                                                  // Initialize max index
            double max_val = vector[i];                                                        // Initialize max value
            for (int j = 1; j < 4; ++j) {                                                    // Find maximum
                if (i + j < vector.size() && vector[i + j] > max_val) {                     // Check if larger
                    max_val = vector[i + j];                                                  // Update max value
                    max_idx = j;                                                              // Update max index
                }
            }
            
            // Decode to base
            char bases[] = {'A', 'T', 'G', 'C'};                                             // Base mapping
            sequence += bases[max_idx];                                                        // Append base
        }
        
        return sequence;                                                                      // Return decoded sequence
    }
    
    double RestrictedBoltzmannMachine::sigmoid(double x) {
        return 1.0 / (1.0 + std::exp(-x));                                                    // Sigmoid function
    }
    
    double RestrictedBoltzmannMachine::sampleBernoulli(double p) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);                                // Uniform distribution
        return dist(rng_) < p ? 1.0 : 0.0;                                                    // Sample from Bernoulli
    }
    
    void RestrictedBoltzmannMachine::contrastiveDivergence(const std::vector<double>& visible, int k) {
        // Positive phase: sample hidden given visible
        std::vector<double> hidden_probs(hidden_size_, 0.0);                                 // Hidden probabilities
        for (int j = 0; j < hidden_size_; ++j) {                                             // Calculate for each hidden unit
            double activation = hidden_bias_[j];                                              // Start with bias
            for (int i = 0; i < visible_size_; ++i) {                                        // Sum over visible
                activation += visible[i] * weights_[i][j];                                   // Weighted sum
            }
            hidden_probs[j] = sigmoid(activation);                                            // Probability
        }
        
        std::vector<double> hidden = hidden_probs;                                            // Use probabilities (mean field)
        
        // Negative phase: k steps of Gibbs sampling
        std::vector<double> visible_neg = visible;                                            // Start from data
        for (int step = 0; step < k; ++step) {                                               // k steps of Gibbs sampling
            // Sample hidden given visible
            for (int j = 0; j < hidden_size_; ++j) {                                         // Sample each hidden unit
                double activation = hidden_bias_[j];                                          // Calculate activation
                for (int i = 0; i < visible_size_; ++i) {                                    // Sum over visible
                    activation += visible_neg[i] * weights_[i][j];                            // Weighted sum
                }
                hidden[j] = sampleBernoulli(sigmoid(activation));                             // Sample from Bernoulli
            }
            
            // Sample visible given hidden
            for (int i = 0; i < visible_size_; ++i) {                                        // Sample each visible unit
                double activation = visible_bias_[i];                                          // Calculate activation
                for (int j = 0; j < hidden_size_; ++j) {                                     // Sum over hidden
                    activation += hidden[j] * weights_[i][j];                                 // Weighted sum
                }
                visible_neg[i] = sampleBernoulli(sigmoid(activation));                        // Sample from Bernoulli
            }
        }
        
        // Update weights and biases
        for (int i = 0; i < visible_size_; ++i) {                                            // Update each weight
            for (int j = 0; j < hidden_size_; ++j) {                                         // Update weight[i][j]
                double positive = visible[i] * hidden_probs[j];                              // Positive phase contribution
                double negative = visible_neg[i] * hidden[j];                                 // Negative phase contribution
                weights_[i][j] += learning_rate_ * (positive - negative);                    // Update weight
            }
            // Update visible bias
            visible_bias_[i] += learning_rate_ * (visible[i] - visible_neg[i]);               // Update visible bias
        }
        
        for (int j = 0; j < hidden_size_; ++j) {                                             // Update hidden biases
            hidden_bias_[j] += learning_rate_ * (hidden_probs[j] - hidden[j]);                // Update hidden bias
        }
    }
    
    // ==================== Deep Belief Network Implementation ====================
    
    DeepBeliefNetwork::DeepBeliefNetwork(const std::vector<int>& layer_sizes)
        : layer_sizes_(layer_sizes) {
        visible_size_ = 64;                                                                    // Default visible size
        
        // Create stack of RBMs
        if (!layer_sizes_.empty()) {                                                          // Check if layers specified
            rbms_.emplace_back(visible_size_, layer_sizes_[0]);                               // First RBM: visible -> first hidden
            for (size_t i = 1; i < layer_sizes_.size(); ++i) {                              // Subsequent RBMs
                rbms_.emplace_back(layer_sizes_[i - 1], layer_sizes_[i]);                    // Previous hidden -> current hidden
            }
        }
    }
    
    void DeepBeliefNetwork::train(const std::vector<std::string>& sequences) {
        if (sequences.empty() || rbms_.empty()) {                                             // Check if valid
            return;                                                                            // Return if invalid
        }
        
        // Greedy layer-wise training
        std::vector<std::string> current_data = sequences;                                   // Start with original sequences
        
        for (size_t layer = 0; layer < rbms_.size(); ++layer) {                              // Train each layer
            rbms_[layer].train(current_data);                                                 // Train current RBM layer
            
            // Encode data through current layer for next layer
            if (layer < rbms_.size() - 1) {                                                   // If not last layer
                std::vector<std::string> encoded_data;                                        // Encoded data for next layer
                for (const auto& seq : current_data) {                                        // Encode each sequence
                    std::vector<double> hidden = rbms_[layer].encode(seq);                     // Get hidden representation
                    // Convert hidden representation back to sequence-like format
                    // (Simplified: in practice would use proper encoding)
                    current_data = sequences;                                                  // Use original for simplicity
                    break;                                                                    // Simplified: use original sequences
                }
            }
        }
    }
    
    std::vector<double> DeepBeliefNetwork::encode(const std::string& sequence) {
        std::vector<double> current = sequenceToVector(sequence);                             // Start with input vector
        
        // Encode through each layer
        for (auto& rbm : rbms_) {                                                             // Pass through each RBM
            // Convert vector to sequence format (simplified: use first RBM's encoding)
            std::string seq;                                                                  // Initialize sequence string
            for (size_t i = 0; i < current.size() && i < 16; ++i) {                          // Convert to sequence (simplified)
                char bases[] = {'A', 'T', 'G', 'C'};
                int idx = static_cast<int>(current[i] * 3.0) % 4;                            // Map to base
                seq += bases[idx];                                                            // Append base
            }
            current = rbm.encode(seq);                                                        // Encode through RBM
        }
        
        return current;                                                                       // Return final representation
    }
    
    std::vector<std::pair<std::string, double>> DeepBeliefNetwork::matchPatterns(
            const std::string& sequence,
            const std::vector<std::string>& pattern_patterns) {
        std::vector<double> seq_features = encode(sequence);                                  // Encode query sequence
        
        std::vector<std::pair<std::string, double>> similarities;                              // Store similarities
        for (const auto& pattern : pattern_patterns) {                                       // Iterate through patterns
            std::vector<double> pat_features = encode(pattern);                               // Encode pattern
            
            // Calculate cosine similarity
            double dot_product = 0.0;                                                         // Initialize dot product
            double norm_seq = 0.0;                                                            // Initialize sequence norm
            double norm_pat = 0.0;                                                            // Initialize pattern norm
            
            for (size_t i = 0; i < seq_features.size() && i < pat_features.size(); ++i) {   // Calculate similarity
                dot_product += seq_features[i] * pat_features[i];                             // Dot product
                norm_seq += seq_features[i] * seq_features[i];                                // Sequence norm squared
                norm_pat += pat_features[i] * pat_features[i];                                 // Pattern norm squared
            }
            
            double similarity = 0.0;                                                          // Initialize similarity
            if (norm_seq > 0.0 && norm_pat > 0.0) {                                          // Avoid division by zero
                similarity = dot_product / (std::sqrt(norm_seq) * std::sqrt(norm_pat));      // Cosine similarity
            }
            
            similarities.push_back({pattern, similarity});                                    // Store similarity
        }
        
        // Sort by similarity
        std::sort(similarities.begin(), similarities.end(),
                  [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
                      return a.second > b.second;                                              // Sort descending
                  });
        
        return similarities;                                                                  // Return sorted similarities
    }
    
    std::vector<double> DeepBeliefNetwork::sequenceToVector(const std::string& sequence) {
        std::vector<double> vector;                                                            // Initialize vector
        vector.reserve(sequence.length() * 4);                                                // Reserve space
        
        for (char c : sequence) {                                                             // Iterate through sequence
            char base = std::toupper(c);                                                      // Convert to uppercase
            vector.push_back(base == 'A' ? 1.0 : 0.0);                                       // One-hot encoding
            vector.push_back(base == 'T' ? 1.0 : 0.0);
            vector.push_back(base == 'G' ? 1.0 : 0.0);
            vector.push_back(base == 'C' ? 1.0 : 0.0);
        }
        
        if (vector.size() < 64) {                                                             // Pad if needed
            vector.resize(64, 0.0);
        } else if (vector.size() > 64) {                                                      // Truncate if needed
            vector.resize(64);
        }
        
        return vector;                                                                        // Return vector
    }
    
    // Helper function for DBN (removed - using inline conversion instead)
    /*
    std::string vectorToSequence(const std::vector<double>& vector) {
        std::string sequence;                                                                 // Initialize sequence
        for (size_t i = 0; i < vector.size(); i += 4) {                                       // Process in groups
            if (i + 3 >= vector.size()) break;                                               // Check bounds
            int max_idx = 0;                                                                  // Find max
            double max_val = vector[i];
            for (int j = 1; j < 4; ++j) {
                if (i + j < vector.size() && vector[i + j] > max_val) {
                    max_val = vector[i + j];
                    max_idx = j;
                }
            }
            char bases[] = {'A', 'T', 'G', 'C'};
            sequence += bases[max_idx];
        }
        return sequence;
    }
    */
    
    // ==================== Markov Decision Process Implementation ====================
    
    MarkovDecisionProcess::MarkovDecisionProcess(double match_reward, double mismatch_penalty,
                                                 double gap_penalty)
        : match_reward_(match_reward), mismatch_penalty_(mismatch_penalty),
          gap_penalty_(gap_penalty), discount_factor_(0.9) {
    }
    
    MarkovDecisionProcess::MDPResult MarkovDecisionProcess::findPattern(
            const std::string& sequence, const std::string& pattern, int max_errors) {
        MDPResult result;                                                                     // Initialize result structure
        
        if (sequence.empty() || pattern.empty()) {                                            // Validate input: check if sequences are non-empty
            return result;                                                                    // Return empty result if invalid input
        }
        
        // Use value iteration to find optimal policy
        std::map<State, double> values = valueIteration(sequence, pattern);                   // Compute optimal value function using value iteration algorithm
        
        // Execute policy: start from initial state
        State current_state;                                                                  // Initialize current state
        current_state.position = 0;                                                           // Start at beginning of sequence (position 0)
        current_state.pattern_index = 0;                                                      // Start at beginning of pattern (index 0)
        current_state.matches = 0;                                                            // Initialize match count to zero
        current_state.errors = max_errors;                                                    // Initialize remaining error budget
        
        double total_reward = 0.0;                                                            // Initialize cumulative reward accumulator
        
        while (current_state.position < static_cast<int>(sequence.length()) &&
               current_state.pattern_index < static_cast<int>(pattern.length())) {            // Continue while valid
            std::vector<Action> actions = getPossibleActions(current_state, sequence, pattern); // Get possible actions
            
            if (actions.empty()) {                                                            // Check if no actions
                break;                                                                        // Break if no actions
            }
            
            // Choose action with maximum Q-value (simplified: greedy policy)
            Action best_action = actions[0];                                                  // Initialize best action
            double best_value = -std::numeric_limits<double>::infinity();                     // Initialize best value
            
            for (const auto& action : actions) {                                              // Find best action
                State next_state = getNextState(current_state, action);                      // Get next state
                double reward = getReward(current_state, action,
                                         sequence[current_state.position],
                                         pattern[current_state.pattern_index]);                // Get reward
                double q_value = reward + discount_factor_ * values[next_state];              // Calculate Q-value
                
                if (q_value > best_value) {                                                   // Check if better
                    best_value = q_value;                                                     // Update best value
                    best_action = action;                                                     // Update best action
                }
            }
            
            // Execute action
            total_reward += getReward(current_state, best_action,
                                     sequence[current_state.position],
                                     pattern[current_state.pattern_index]);                    // Add reward
            result.action_sequence.push_back(best_action);                                    // Record action
            current_state = getNextState(current_state, best_action);                         // Update state
            
            if (best_action.type == Action::MATCH && current_state.matches > 0) {            // Check if match found
                result.positions.push_back(current_state.position - 1);                      // Record position
            }
        }
        
        result.total_reward = total_reward;                                                   // Store total reward
        return result;                                                                        // Return result
    }
    
    std::map<MarkovDecisionProcess::State, double> MarkovDecisionProcess::valueIteration(
            const std::string& sequence, const std::string& pattern) {
        std::map<State, double> values;                                                       // Initialize value function
        const double epsilon = 0.001;                                                         // Convergence threshold
        const int max_iterations = 1000;                                                      // Maximum iterations
        
        // Initialize values
        for (int i = 0; i <= static_cast<int>(sequence.length()); ++i) {                     // Initialize all states
            for (int j = 0; j <= static_cast<int>(pattern.length()); ++j) {
                State state;
                state.position = i;
                state.pattern_index = j;
                state.matches = 0;
                state.errors = 0;
                values[state] = 0.0;                                                          // Initialize to zero
            }
        }
        
        // Value iteration
        for (int iter = 0; iter < max_iterations; ++iter) {                                  // Iterate until convergence
            std::map<State, double> new_values = values;                                      // Copy current values
            double max_change = 0.0;                                                          // Track maximum change
            
            for (auto& value_pair : values) {                                                 // Update each state
                State state = value_pair.first;                                               // Get state
                if (state.position >= static_cast<int>(sequence.length()) ||
                    state.pattern_index >= static_cast<int>(pattern.length())) {              // Terminal state
                    continue;                                                                 // Skip terminal states
                }
                
                std::vector<Action> actions = getPossibleActions(state, sequence, pattern);    // Get possible actions
                double max_q = -std::numeric_limits<double>::infinity();                     // Initialize max Q-value
                
                for (const auto& action : actions) {                                          // Find maximum Q-value
                    State next_state = getNextState(state, action);                          // Get next state
                    double reward = getReward(state, action,
                                             sequence[state.position],
                                             pattern[state.pattern_index]);                    // Get reward
                    double q_value = reward + discount_factor_ * values[next_state];         // Calculate Q-value
                    max_q = std::max(max_q, q_value);                                         // Update maximum
                }
                
                new_values[state] = max_q;                                                     // Update value
                max_change = std::max(max_change, std::abs(new_values[state] - values[state])); // Track change
            }
            
            values = new_values;                                                               // Update values
            if (max_change < epsilon) {                                                       // Check convergence
                break;                                                                        // Converged
            }
        }
        
        return values;                                                                        // Return value function
    }
    
    double MarkovDecisionProcess::getReward(const State& state, const Action& action,
                                            char seq_char, char pat_char) {
        switch (action.type) {                                                                // Switch on action type
            case Action::MATCH:                                                               // Match action
                if (std::toupper(seq_char) == std::toupper(pat_char)) {                      // Check if characters match
                    return match_reward_;                                                     // Return match reward
                } else {
                    return mismatch_penalty_;                                                  // Return mismatch penalty
                }
            case Action::INSERT:                                                              // Insert action
            case Action::DELETE:                                                              // Delete action
                return gap_penalty_;                                                          // Return gap penalty
            case Action::SKIP:                                                                // Skip action
                return -0.5;                                                                  // Small penalty for skipping
            default:
                return 0.0;                                                                   // Default reward
        }
    }
    
    MarkovDecisionProcess::State MarkovDecisionProcess::getNextState(const State& current, const Action& action) {
        State next = current;                                                                 // Copy current state
        
        switch (action.type) {                                                                // Switch on action type
            case Action::MATCH:                                                               // Match: advance both
                next.position++;
                next.pattern_index++;
                next.matches++;
                break;
            case Action::INSERT:                                                              // Insert: advance sequence only
                next.position++;
                break;
            case Action::DELETE:                                                              // Delete: advance pattern only
                next.pattern_index++;
                next.errors--;
                break;
            case Action::SKIP:                                                                // Skip: advance sequence only
                next.position++;
                break;
        }
        
        return next;                                                                          // Return next state
    }
    
    std::vector<MarkovDecisionProcess::Action> MarkovDecisionProcess::getPossibleActions(
            const State& state, const std::string& sequence, const std::string& pattern) {
        std::vector<Action> actions;                                                          // Initialize actions
        
        if (state.position < static_cast<int>(sequence.length()) &&
            state.pattern_index < static_cast<int>(pattern.length())) {                       // If both valid
            actions.push_back({Action::MATCH, 0.0});                                         // Can match
            actions.push_back({Action::INSERT, 0.0});                                        // Can insert
            actions.push_back({Action::DELETE, 0.0});                                        // Can delete
        } else if (state.position < static_cast<int>(sequence.length())) {                    // If only sequence valid
            actions.push_back({Action::SKIP, 0.0});                                          // Can skip
        } else if (state.pattern_index < static_cast<int>(pattern.length())) {                // If only pattern valid
            actions.push_back({Action::DELETE, 0.0});                                        // Can delete
        }
        
        return actions;                                                                       // Return possible actions
    }
    
    // ==================== Markov Random Field Implementation ====================
    
    MarkovRandomField::MarkovRandomField(double match_energy, double mismatch_energy,
                                        double smoothness_weight)
        : match_energy_(match_energy), mismatch_energy_(mismatch_energy),
          smoothness_weight_(smoothness_weight), rng_(std::random_device{}()) {
    }
    
    MarkovRandomField::MRFResult MarkovRandomField::findPattern(
            const std::string& sequence, const std::string& pattern) {
        return gibbsSampling(sequence, pattern, 1000);                                        // Use Gibbs sampling for inference
    }
    
    MarkovRandomField::MRFResult MarkovRandomField::gibbsSampling(
            const std::string& sequence, const std::string& pattern, int num_iterations) {
        MRFResult result;                                                                     // Initialize result
        
        if (sequence.length() < pattern.length()) {                                           // Check if valid
            return result;                                                                    // Return empty if invalid
        }
        
        int num_positions = static_cast<int>(sequence.length() - pattern.length() + 1);       // Number of possible positions
        std::vector<bool> matches(num_positions, false);                                      // Initialize match variables
        
        // Initialize randomly
        std::uniform_real_distribution<double> init_dist(0.0, 1.0);                          // Random initialization
        for (int i = 0; i < num_positions; ++i) {                                            // Initialize each position
            matches[i] = init_dist(rng_) > 0.5;                                              // Random initialization
        }
        
        // Gibbs sampling iterations
        for (int iter = 0; iter < num_iterations; ++iter) {                                  // Iterate for burn-in and sampling
            for (int i = 0; i < num_positions; ++i) {                                        // Update each variable
                matches[i] = sampleConditional(i, matches, sequence, pattern);                // Sample from conditional
            }
        }
        
        // Collect results
        for (int i = 0; i < num_positions; ++i) {                                             // Record matches
            if (matches[i]) {                                                                 // If position matches
                result.matched_positions.push_back(i);                                       // Record position
            }
        }
        
        result.energy = calculateEnergy(matches, sequence, pattern);                         // Calculate energy
        
        return result;                                                                        // Return result
    }
    
    MarkovRandomField::MRFResult MarkovRandomField::beliefPropagation(
            const std::string& sequence, const std::string& pattern) {
        MRFResult result;                                                                     // Initialize result
        
        if (sequence.length() < pattern.length()) {                                           // Check if valid
            return result;                                                                    // Return empty if invalid
        }
        
        int num_positions = static_cast<int>(sequence.length() - pattern.length() + 1);       // Number of positions
        
        // Initialize messages (simplified belief propagation)
        std::map<int, double> beliefs;                                                        // Initialize beliefs
        
        for (int i = 0; i < num_positions; ++i) {                                            // Calculate belief for each position
            // Calculate local evidence
            double local_energy = 0.0;                                                        // Initialize local energy
            for (size_t j = 0; j < pattern.length() && (i + j) < sequence.length(); ++j) {   // Check match quality
                if (std::toupper(sequence[i + j]) == std::toupper(pattern[j])) {            // If characters match
                    local_energy += match_energy_;                                             // Add match energy
                } else {
                    local_energy += mismatch_energy_;                                         // Add mismatch energy
                }
            }
            
            // Convert energy to probability (Boltzmann distribution)
            double prob = 1.0 / (1.0 + std::exp(local_energy));                               // Sigmoid of negative energy
            beliefs[i] = prob;                                                                // Store belief
            
            if (prob > 0.5) {                                                                 // If probability > threshold
                result.matched_positions.push_back(i);                                       // Record as match
                result.position_probabilities[i] = prob;                                       // Store probability
            }
        }
        
        result.energy = 0.0;                                                                  // Calculate total energy (simplified)
        for (const auto& belief_pair : beliefs) {                                             // Sum energies
            if (belief_pair.second > 0.5) {                                                   // If matched
                result.energy += mismatch_energy_;                                            // Add energy (simplified)
            }
        }
        
        return result;                                                                        // Return result
    }
    
    double MarkovRandomField::calculateEnergy(const std::vector<bool>& matches,
                                              const std::string& sequence, const std::string& pattern) {
        double energy = 0.0;                                                                  // Initialize energy
        
        for (size_t i = 0; i < matches.size(); ++i) {                                        // Iterate through positions
            if (matches[i]) {                                                                 // If position matches
                // Local energy: match quality
                for (size_t j = 0; j < pattern.length() && (i + j) < sequence.length(); ++j) { // Check each character
                    if (std::toupper(sequence[i + j]) == std::toupper(pattern[j])) {        // If characters match
                        energy += match_energy_;                                               // Add match energy
                    } else {
                        energy += mismatch_energy_;                                           // Add mismatch energy
                    }
                }
                
                // Smoothness energy: encourage nearby matches
                if (i > 0 && matches[i - 1]) {                                               // If previous position matches
                    energy -= smoothness_weight_;                                             // Reduce energy (encourage smoothness)
                }
                if (i + 1 < matches.size() && matches[i + 1]) {                              // If next position matches
                    energy -= smoothness_weight_;                                             // Reduce energy
                }
            }
        }
        
        return energy;                                                                        // Return total energy
    }
    
    double MarkovRandomField::localEnergy(int position, bool is_match,
                                         const std::vector<bool>& neighbors,
                                         const std::string& sequence, const std::string& pattern) {
        double energy = 0.0;                                                                  // Initialize energy
        
        if (is_match) {                                                                       // If this position matches
            // Match quality energy
            for (size_t j = 0; j < pattern.length() && (position + j) < sequence.length(); ++j) { // Check match quality
                if (std::toupper(sequence[position + j]) == std::toupper(pattern[j])) {     // If characters match
                    energy += match_energy_;                                                   // Add match energy
                } else {
                    energy += mismatch_energy_;                                               // Add mismatch energy
                }
            }
            
            // Neighbor smoothness
            for (bool neighbor : neighbors) {                                                  // Check neighbors
                if (neighbor) {                                                               // If neighbor matches
                    energy -= smoothness_weight_;                                             // Reduce energy (smoothness)
                }
            }
        }
        
        return energy;                                                                        // Return local energy
    }
    
    bool MarkovRandomField::sampleConditional(int position, const std::vector<bool>& current_state,
                                              const std::string& sequence, const std::string& pattern) {
        // Get neighbors
        std::vector<bool> neighbors;                                                          // Initialize neighbors
        if (position > 0) {                                                                  // If not first position
            neighbors.push_back(current_state[position - 1]);                                 // Add previous neighbor
        }
        if (position + 1 < static_cast<int>(current_state.size())) {                         // If not last position
            neighbors.push_back(current_state[position + 1]);                                 // Add next neighbor
        }
        
        // Calculate energy for both states
        double energy_match = localEnergy(position, true, neighbors, sequence, pattern);      // Energy if match
        double energy_no_match = localEnergy(position, false, neighbors, sequence, pattern); // Energy if no match
        
        // Convert to probabilities (Boltzmann distribution)
        double prob_match = 1.0 / (1.0 + std::exp(energy_match - energy_no_match));          // Probability of match
        
        // Sample from Bernoulli
        std::uniform_real_distribution<double> dist(0.0, 1.0);                               // Uniform distribution
        return dist(rng_) < prob_match;                                                       // Sample and return
    }
    
} // namespace ProbabilisticML

