#include "cnn_sequence.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <cctype>

CNNSequenceModel::CNNSequenceModel(int input_size, int num_filters, int filter_size)
    : input_size_(input_size), num_filters_(num_filters), filter_size_(filter_size) {
    initializeWeights();
}

double CNNSequenceModel::predict(const std::vector<double>& sequence) {
    if (sequence.size() != static_cast<size_t>(input_size_)) {
        return 0.0;
    }
    
    // Extract features
    std::vector<double> features = extractFeatures(sequence);
    
    // Dense layer
    double output = 0.0;
    for (size_t i = 0; i < features.size() && i < dense_weights_[0].size(); ++i) {
        output += features[i] * dense_weights_[0][i];
    }
    output += dense_bias_[0];
    
    // Sigmoid activation
    return sigmoid(output);
}

std::vector<double> CNNSequenceModel::encodeSequence(const std::string& sequence) {
    std::vector<double> encoded;
    
    for (char c : sequence) {
        char base = std::toupper(c);
        if (base == 'A') encoded.push_back(0.0);
        else if (base == 'T') encoded.push_back(1.0);
        else if (base == 'G') encoded.push_back(2.0);
        else if (base == 'C') encoded.push_back(3.0);
        else encoded.push_back(0.0);  // Unknown base
    }
    
    // Normalize to [0, 1]
    for (double& val : encoded) {
        val /= 3.0;
    }
    
    // Pad or truncate to fixed size (for static function)
    const int fixed_size = 128;  // Default encoding size
    if (encoded.size() < static_cast<size_t>(fixed_size)) {
        encoded.resize(fixed_size, 0.0);
    } else if (encoded.size() > static_cast<size_t>(fixed_size)) {
        encoded.resize(fixed_size);
    }
    
    return encoded;
}

std::vector<double> CNNSequenceModel::extractFeatures(const std::vector<double>& sequence) {
    std::vector<double> features;
    
    // Apply convolutional filters
    for (int f = 0; f < num_filters_; ++f) {
        std::vector<double> conv_output = convolve(sequence, conv_weights_[f]);
        
        // Apply ReLU
        for (double& val : conv_output) {
            val = relu(val);
        }
        
        // Max pooling
        std::vector<double> pooled = maxPooling(conv_output, 2);
        
        // Flatten and add to features
        for (double val : pooled) {
            features.push_back(val);
        }
    }
    
    return features;
}

std::vector<size_t> CNNSequenceModel::searchPattern(const std::string& sequence,
                                                    const std::string& pattern,
                                                    double threshold) {
    std::vector<size_t> positions;
    
    if (pattern.length() > sequence.length()) {
        return positions;
    }
    
    // Slide window over sequence
    for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {
        std::string window = sequence.substr(i, pattern.length());
        
        // Encode window
        std::vector<double> encoded = encodeSequence(window);
        
        // Predict match probability
        double probability = predict(encoded);
        
        if (probability >= threshold) {
            positions.push_back(i);
        }
    }
    
    return positions;
}

void CNNSequenceModel::train(const std::vector<std::string>& sequences,
                            const std::vector<double>& labels) {
    // Simplified training - in real implementation, would use backpropagation
    // This is a placeholder that could be extended with actual gradient descent
    
    if (sequences.size() != labels.size()) {
        return;
    }
    
    // For demonstration, we'll do a simple update
    // In practice, this would involve:
    // 1. Forward pass
    // 2. Calculate loss
    // 3. Backpropagation
    // 4. Weight updates
    
    // Placeholder: weights are initialized randomly
    // Real training would update them based on gradients
}

std::vector<double> CNNSequenceModel::convolve(const std::vector<double>& input,
                                              const std::vector<double>& filter) {
    std::vector<double> output;
    
    if (input.size() < filter.size()) {
        return output;
    }
    
    for (size_t i = 0; i <= input.size() - filter.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < filter.size(); ++j) {
            sum += input[i + j] * filter[j];
        }
        output.push_back(sum);
    }
    
    return output;
}

double CNNSequenceModel::relu(double x) {
    return std::max(0.0, x);
}

double CNNSequenceModel::sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

std::vector<double> CNNSequenceModel::maxPooling(const std::vector<double>& input, int pool_size) {
    std::vector<double> output;
    
    for (size_t i = 0; i < input.size(); i += pool_size) {
        double max_val = input[i];
        for (int j = 1; j < pool_size && (i + j) < input.size(); ++j) {
            max_val = std::max(max_val, input[i + j]);
        }
        output.push_back(max_val);
    }
    
    return output;
}

void CNNSequenceModel::initializeWeights() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-0.1, 0.1);
    
    // Initialize convolutional filters
    conv_weights_.resize(num_filters_);
    conv_bias_.resize(num_filters_);
    
    for (int i = 0; i < num_filters_; ++i) {
        conv_weights_[i].resize(filter_size_);
        for (int j = 0; j < filter_size_; ++j) {
            conv_weights_[i][j] = dist(gen);
        }
        conv_bias_[i] = dist(gen);
    }
    
    // Initialize dense layer (simplified - single output)
    int feature_size = (input_size_ - filter_size_ + 1) / 2 * num_filters_;
    dense_weights_.resize(1);
    dense_weights_[0].resize(feature_size);
    for (size_t i = 0; i < dense_weights_[0].size(); ++i) {
        dense_weights_[0][i] = dist(gen);
    }
    dense_bias_.resize(1);
    dense_bias_[0] = dist(gen);
}

