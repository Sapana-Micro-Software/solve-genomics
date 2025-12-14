#ifndef DEEP_LEARNING_ADVANCED_H
#define DEEP_LEARNING_ADVANCED_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <memory>

/**
 * Advanced Deep Learning Approaches for DNA Sequence Analysis
 */
class DeepLearningAdvanced {
public:
    /**
     * LSTM (Long Short-Term Memory) for sequence modeling
     */
    class LSTM {
    public:
        LSTM(int input_size, int hidden_size, int num_layers = 1);
        
        /**
         * Process sequence through LSTM
         * @param sequence Input sequence
         * @return Hidden state vector
         */
        std::vector<double> forward(const std::string& sequence);
        
        /**
         * Search pattern using LSTM embeddings
         * @param sequence Sequence to search
         * @param pattern Pattern to find
         * @param threshold Similarity threshold
         * @return Positions and similarities
         */
        std::vector<std::pair<size_t, double>> searchPattern(const std::string& sequence,
                                                            const std::string& pattern,
                                                            double threshold = 0.5);
        
    private:
        int input_size_;
        int hidden_size_;
        int num_layers_;
        std::vector<std::vector<double>> weights_;
        std::vector<double> biases_;
        
        double sigmoid(double x);
        double tanh_activation(double x);
        std::vector<double> encodeSequence(const std::string& sequence);
    };
    
    /**
     * GRU (Gated Recurrent Unit) for sequence processing
     */
    class GRU {
    public:
        GRU(int input_size, int hidden_size);
        
        std::vector<double> forward(const std::string& sequence);
        std::vector<std::pair<size_t, double>> searchPattern(const std::string& sequence,
                                                            const std::string& pattern,
                                                            double threshold = 0.5);
        
    private:
        int input_size_;
        int hidden_size_;
        std::vector<std::vector<double>> weights_;
        std::vector<double> biases_;
        
        double sigmoid(double x);
        double tanh_activation(double x);
        std::vector<double> encodeSequence(const std::string& sequence);
    };
    
    /**
     * Attention mechanism for sequence alignment
     */
    class AttentionModel {
    public:
        AttentionModel(int embedding_dim, int num_heads = 4);
        
        /**
         * Compute attention-weighted sequence representation
         * @param sequence Input sequence
         * @return Attention-weighted embedding
         */
        std::vector<double> computeAttention(const std::string& sequence);
        
        /**
         * Align two sequences using attention
         * @param seq1 First sequence
         * @param seq2 Second sequence
         * @return Alignment score
         */
        double alignSequences(const std::string& seq1, const std::string& seq2);
        
    private:
        int embedding_dim_;
        int num_heads_;
        std::vector<std::vector<double>> query_weights_;
        std::vector<std::vector<double>> key_weights_;
        std::vector<std::vector<double>> value_weights_;
        
        std::vector<double> encodeSequence(const std::string& sequence);
        std::vector<std::vector<double>> multiHeadAttention(const std::vector<std::vector<double>>& queries,
                                                            const std::vector<std::vector<double>>& keys,
                                                            const std::vector<std::vector<double>>& values);
    };
    
    /**
     * Siamese Network for sequence similarity
     */
    class SiameseNetwork {
    public:
        SiameseNetwork(int embedding_dim);
        
        /**
         * Compute similarity between two sequences
         * @param seq1 First sequence
         * @param seq2 Second sequence
         * @return Similarity score (0-1)
         */
        double computeSimilarity(const std::string& seq1, const std::string& seq2);
        
        /**
         * Find similar sequences
         * @param query Query sequence
         * @param sequences Database sequences
         * @param top_k Number of results
         * @return Top-k similar sequences with scores
         */
        std::vector<std::pair<size_t, double>> findSimilar(const std::string& query,
                                                          const std::vector<std::string>& sequences,
                                                          int top_k = 10);
        
    private:
        int embedding_dim_;
        std::vector<std::vector<double>> shared_weights_;
        
        std::vector<double> encodeSequence(const std::string& sequence);
        std::vector<double> forwardPass(const std::string& sequence);
    };
};

#endif // DEEP_LEARNING_ADVANCED_H

