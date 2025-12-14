#include "deep_learning_advanced.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <cctype>

// LSTM Implementation

DeepLearningAdvanced::LSTM::LSTM(int input_size, int hidden_size, int num_layers)
    : input_size_(input_size), hidden_size_(hidden_size), num_layers_(num_layers) {
    // Initialize weights randomly (simplified)
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 0.1);
    
    int total_params = (input_size_ + hidden_size_ + 1) * hidden_size_ * 4;  // 4 gates
    weights_.resize(num_layers_);
    biases_.resize(hidden_size_ * 4);
    
    for (int layer = 0; layer < num_layers_; ++layer) {
        weights_[layer].resize(total_params);
        for (double& w : weights_[layer]) {
            w = dist(rng);
        }
    }
    
    for (double& b : biases_) {
        b = dist(rng);
    }
}

double DeepLearningAdvanced::LSTM::sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

double DeepLearningAdvanced::LSTM::tanh_activation(double x) {
    return std::tanh(x);
}

std::vector<double> DeepLearningAdvanced::LSTM::encodeSequence(const std::string& sequence) {
    std::vector<double> encoded;
    encoded.reserve(sequence.length() * input_size_);
    
    for (char c : sequence) {
        char upper_c = std::toupper(c);
        std::vector<double> char_vec(input_size_, 0.0);
        
        // One-hot encoding simplified
        if (upper_c == 'A') char_vec[0] = 1.0;
        else if (upper_c == 'T') char_vec[1] = 1.0;
        else if (upper_c == 'G') char_vec[2] = 1.0;
        else if (upper_c == 'C') char_vec[3] = 1.0;
        
        encoded.insert(encoded.end(), char_vec.begin(), char_vec.end());
    }
    
    return encoded;
}

std::vector<double> DeepLearningAdvanced::LSTM::forward(const std::string& sequence) {
    std::vector<double> encoded = encodeSequence(sequence);
    std::vector<double> hidden(hidden_size_, 0.0);
    std::vector<double> cell(hidden_size_, 0.0);
    
    // Simplified LSTM forward pass
    for (size_t t = 0; t < sequence.length(); ++t) {
        // Extract input at time t
        std::vector<double> input(input_size_);
        for (int i = 0; i < input_size_; ++i) {
            input[i] = encoded[t * input_size_ + i];
        }
        
        // LSTM gates (simplified computation)
        std::vector<double> forget_gate(hidden_size_, 0.5);  // Simplified
        std::vector<double> input_gate(hidden_size_, 0.5);
        std::vector<double> output_gate(hidden_size_, 0.5);
        std::vector<double> candidate(hidden_size_, 0.0);
        
        // Update cell and hidden state
        for (int i = 0; i < hidden_size_; ++i) {
            cell[i] = forget_gate[i] * cell[i] + input_gate[i] * candidate[i];
            hidden[i] = output_gate[i] * tanh_activation(cell[i]);
        }
    }
    
    return hidden;
}

std::vector<std::pair<size_t, double>> DeepLearningAdvanced::LSTM::searchPattern(
    const std::string& sequence, const std::string& pattern, double threshold) {
    
    std::vector<std::pair<size_t, double>> results;
    
    std::vector<double> pattern_embedding = forward(pattern);
    
    for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {
        std::string window = sequence.substr(i, pattern.length());
        std::vector<double> window_embedding = forward(window);
        
        // Cosine similarity
        double dot = 0.0, norm1 = 0.0, norm2 = 0.0;
        for (size_t j = 0; j < pattern_embedding.size(); ++j) {
            dot += pattern_embedding[j] * window_embedding[j];
            norm1 += pattern_embedding[j] * pattern_embedding[j];
            norm2 += window_embedding[j] * window_embedding[j];
        }
        
        double similarity = dot / (std::sqrt(norm1) * std::sqrt(norm2) + 1e-8);
        
        if (similarity >= threshold) {
            results.push_back({i, similarity});
        }
    }
    
    return results;
}

// GRU Implementation

DeepLearningAdvanced::GRU::GRU(int input_size, int hidden_size)
    : input_size_(input_size), hidden_size_(hidden_size) {
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 0.1);
    
    int total_params = (input_size_ + hidden_size_ + 1) * hidden_size_ * 3;  // 3 gates
    weights_.resize(1);
    weights_[0].resize(total_params);
    biases_.resize(hidden_size_ * 3);
    
    for (double& w : weights_[0]) {
        w = dist(rng);
    }
    for (double& b : biases_) {
        b = dist(rng);
    }
}

double DeepLearningAdvanced::GRU::sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

double DeepLearningAdvanced::GRU::tanh_activation(double x) {
    return std::tanh(x);
}

std::vector<double> DeepLearningAdvanced::GRU::encodeSequence(const std::string& sequence) {
    std::vector<double> encoded;
    for (char c : sequence) {
        char upper_c = std::toupper(c);
        if (upper_c == 'A') encoded.push_back(1.0);
        else if (upper_c == 'T') encoded.push_back(0.0);
        else if (upper_c == 'G') encoded.push_back(0.5);
        else if (upper_c == 'C') encoded.push_back(0.25);
        else encoded.push_back(0.0);
    }
    return encoded;
}

std::vector<double> DeepLearningAdvanced::GRU::forward(const std::string& sequence) {
    std::vector<double> encoded = encodeSequence(sequence);
    std::vector<double> hidden(hidden_size_, 0.0);
    
    for (size_t t = 0; t < sequence.length(); ++t) {
        // Simplified GRU computation
        std::vector<double> reset_gate(hidden_size_, 0.5);
        std::vector<double> update_gate(hidden_size_, 0.5);
        std::vector<double> candidate(hidden_size_, 0.0);
        
        for (int i = 0; i < hidden_size_; ++i) {
            hidden[i] = (1.0 - update_gate[i]) * hidden[i] + update_gate[i] * candidate[i];
        }
    }
    
    return hidden;
}

std::vector<std::pair<size_t, double>> DeepLearningAdvanced::GRU::searchPattern(
    const std::string& sequence, const std::string& pattern, double threshold) {
    
    std::vector<std::pair<size_t, double>> results;
    std::vector<double> pattern_embedding = forward(pattern);
    
    for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {
        std::string window = sequence.substr(i, pattern.length());
        std::vector<double> window_embedding = forward(window);
        
        double dot = 0.0, norm1 = 0.0, norm2 = 0.0;
        for (size_t j = 0; j < pattern_embedding.size(); ++j) {
            dot += pattern_embedding[j] * window_embedding[j];
            norm1 += pattern_embedding[j] * pattern_embedding[j];
            norm2 += window_embedding[j] * window_embedding[j];
        }
        
        double similarity = dot / (std::sqrt(norm1) * std::sqrt(norm2) + 1e-8);
        
        if (similarity >= threshold) {
            results.push_back({i, similarity});
        }
    }
    
    return results;
}

// Attention Model Implementation

DeepLearningAdvanced::AttentionModel::AttentionModel(int embedding_dim, int num_heads)
    : embedding_dim_(embedding_dim), num_heads_(num_heads) {
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 0.1);
    
    query_weights_.resize(num_heads_);
    key_weights_.resize(num_heads_);
    value_weights_.resize(num_heads_);
    
    for (int h = 0; h < num_heads_; ++h) {
        query_weights_[h].resize(embedding_dim_ * embedding_dim_);
        key_weights_[h].resize(embedding_dim_ * embedding_dim_);
        value_weights_[h].resize(embedding_dim_ * embedding_dim_);
        
        for (double& w : query_weights_[h]) w = dist(rng);
        for (double& w : key_weights_[h]) w = dist(rng);
        for (double& w : value_weights_[h]) w = dist(rng);
    }
}

std::vector<double> DeepLearningAdvanced::AttentionModel::encodeSequence(const std::string& sequence) {
    std::vector<double> embedding(embedding_dim_, 0.0);
    
    // Simple encoding: frequency-based
    int counts[4] = {0, 0, 0, 0};  // A, T, G, C
    for (char c : sequence) {
        char upper_c = std::toupper(c);
        if (upper_c == 'A') counts[0]++;
        else if (upper_c == 'T') counts[1]++;
        else if (upper_c == 'G') counts[2]++;
        else if (upper_c == 'C') counts[3]++;
    }
    
    int total = static_cast<int>(sequence.length());
    if (total > 0) {
        for (int i = 0; i < 4 && i < embedding_dim_; ++i) {
            embedding[i] = static_cast<double>(counts[i]) / total;
        }
    }
    
    return embedding;
}

std::vector<std::vector<double>> DeepLearningAdvanced::AttentionModel::multiHeadAttention(
    const std::vector<std::vector<double>>& queries,
    const std::vector<std::vector<double>>& keys,
    const std::vector<std::vector<double>>& values) {
    
    std::vector<std::vector<double>> output(queries.size(), std::vector<double>(embedding_dim_, 0.0));
    
    // Simplified multi-head attention
    for (size_t i = 0; i < queries.size(); ++i) {
        for (int h = 0; h < num_heads_; ++h) {
            // Compute attention scores
            double score = 0.0;
            for (size_t j = 0; j < keys.size(); ++j) {
                double dot = 0.0;
                for (size_t k = 0; k < queries[i].size() && k < keys[j].size(); ++k) {
                    dot += queries[i][k] * keys[j][k];
                }
                score += std::exp(dot);
            }
            
            // Weighted sum
            for (size_t j = 0; j < values.size(); ++j) {
                double weight = std::exp(score) / (score + 1e-8);
                for (size_t k = 0; k < output[i].size() && k < values[j].size(); ++k) {
                    output[i][k] += weight * values[j][k];
                }
            }
        }
    }
    
    return output;
}

std::vector<double> DeepLearningAdvanced::AttentionModel::computeAttention(const std::string& sequence) {
    std::vector<double> embedding = encodeSequence(sequence);
    
    // Create queries, keys, values from embedding
    std::vector<std::vector<double>> queries = {embedding};
    std::vector<std::vector<double>> keys = {embedding};
    std::vector<std::vector<double>> values = {embedding};
    
    std::vector<std::vector<double>> attended = multiHeadAttention(queries, keys, values);
    
    return attended.empty() ? embedding : attended[0];
}

double DeepLearningAdvanced::AttentionModel::alignSequences(const std::string& seq1, const std::string& seq2) {
    std::vector<double> emb1 = computeAttention(seq1);
    std::vector<double> emb2 = computeAttention(seq2);
    
    double dot = 0.0, norm1 = 0.0, norm2 = 0.0;
    for (size_t i = 0; i < emb1.size() && i < emb2.size(); ++i) {
        dot += emb1[i] * emb2[i];
        norm1 += emb1[i] * emb1[i];
        norm2 += emb2[i] * emb2[i];
    }
    
    return dot / (std::sqrt(norm1) * std::sqrt(norm2) + 1e-8);
}

// Siamese Network Implementation

DeepLearningAdvanced::SiameseNetwork::SiameseNetwork(int embedding_dim) : embedding_dim_(embedding_dim) {
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 0.1);
    
    shared_weights_.resize(embedding_dim_);
    for (auto& row : shared_weights_) {
        row.resize(embedding_dim_);
        for (double& w : row) {
            w = dist(rng);
        }
    }
}

std::vector<double> DeepLearningAdvanced::SiameseNetwork::encodeSequence(const std::string& sequence) {
    std::vector<double> embedding(embedding_dim_, 0.0);
    
    // K-mer based encoding
    int k = 3;
    for (size_t i = 0; i <= sequence.length() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        // Simple hash-based embedding
        size_t hash = std::hash<std::string>{}(kmer);
        embedding[hash % embedding_dim_] += 1.0;
    }
    
    // Normalize
    double norm = 0.0;
    for (double val : embedding) norm += val * val;
    norm = std::sqrt(norm);
    if (norm > 0) {
        for (double& val : embedding) val /= norm;
    }
    
    return embedding;
}

std::vector<double> DeepLearningAdvanced::SiameseNetwork::forwardPass(const std::string& sequence) {
    std::vector<double> encoded = encodeSequence(sequence);
    std::vector<double> output(embedding_dim_, 0.0);
    
    // Shared network forward pass
    for (int i = 0; i < embedding_dim_; ++i) {
        for (int j = 0; j < embedding_dim_; ++j) {
            output[i] += shared_weights_[i][j] * encoded[j];
        }
        output[i] = std::tanh(output[i]);  // Activation
    }
    
    return output;
}

double DeepLearningAdvanced::SiameseNetwork::computeSimilarity(const std::string& seq1, const std::string& seq2) {
    std::vector<double> emb1 = forwardPass(seq1);
    std::vector<double> emb2 = forwardPass(seq2);
    
    double dot = 0.0, norm1 = 0.0, norm2 = 0.0;
    for (size_t i = 0; i < emb1.size() && i < emb2.size(); ++i) {
        dot += emb1[i] * emb2[i];
        norm1 += emb1[i] * emb1[i];
        norm2 += emb2[i] * emb2[i];
    }
    
    return dot / (std::sqrt(norm1) * std::sqrt(norm2) + 1e-8);
}

std::vector<std::pair<size_t, double>> DeepLearningAdvanced::SiameseNetwork::findSimilar(
    const std::string& query, const std::vector<std::string>& sequences, int top_k) {
    
    std::vector<std::pair<size_t, double>> similarities;
    
    for (size_t i = 0; i < sequences.size(); ++i) {
        double sim = computeSimilarity(query, sequences[i]);
        similarities.push_back({i, sim});
    }
    
    // Sort by similarity
    std::sort(similarities.begin(), similarities.end(),
        [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
            return a.second > b.second;
        });
    
    if (static_cast<int>(similarities.size()) > top_k) {
        similarities.resize(top_k);
    }
    
    return similarities;
}

