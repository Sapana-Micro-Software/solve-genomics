#ifndef CNN_SEQUENCE_H
#define CNN_SEQUENCE_H

#include <string>
#include <vector>
#include <memory>

/**
 * Lightweight CNN for DNA sequence analysis
 * Simplified convolutional neural network for pattern recognition
 */
class CNNSequenceModel {
public:
    CNNSequenceModel(int input_size = 100, int num_filters = 32, int filter_size = 3);
    
    /**
     * Forward pass: predict pattern match probability
     * @param sequence Input DNA sequence (encoded)
     * @return Probability of pattern match (0.0 to 1.0)
     */
    double predict(const std::vector<double>& sequence);
    
    /**
     * Encode DNA sequence to numeric vector
     * @param sequence DNA sequence string
     * @return Encoded vector
     */
    static std::vector<double> encodeSequence(const std::string& sequence);
    
    /**
     * Extract features using convolutional layers
     * @param sequence Encoded sequence
     * @return Feature vector
     */
    std::vector<double> extractFeatures(const std::vector<double>& sequence);
    
    /**
     * Search for pattern using CNN
     * @param sequence DNA sequence to search in
     * @param pattern Pattern to find
     * @param threshold Probability threshold
     * @return Positions where pattern likely matches
     */
    std::vector<size_t> searchPattern(const std::string& sequence, 
                                     const std::string& pattern,
                                     double threshold = 0.5);
    
    /**
     * Train model (simplified - would need actual training data)
     * @param sequences Training sequences
     * @param labels Training labels (1.0 for match, 0.0 for no match)
     */
    void train(const std::vector<std::string>& sequences, 
              const std::vector<double>& labels);
    
private:
    int input_size_;
    int num_filters_;
    int filter_size_;
    
    // Simplified CNN weights (in real implementation, these would be learned)
    std::vector<std::vector<double>> conv_weights_;
    std::vector<double> conv_bias_;
    std::vector<std::vector<double>> dense_weights_;
    std::vector<double> dense_bias_;
    
    /**
     * Convolution operation
     */
    std::vector<double> convolve(const std::vector<double>& input, 
                                const std::vector<double>& filter);
    
    /**
     * ReLU activation
     */
    double relu(double x);
    
    /**
     * Sigmoid activation
     */
    double sigmoid(double x);
    
    /**
     * Max pooling
     */
    std::vector<double> maxPooling(const std::vector<double>& input, int pool_size);
    
    /**
     * Initialize weights (random initialization)
     */
    void initializeWeights();
};

#endif // CNN_SEQUENCE_H

