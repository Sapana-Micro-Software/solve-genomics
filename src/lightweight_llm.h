#ifndef LIGHTWEIGHT_LLM_H
#define LIGHTWEIGHT_LLM_H

#include <string>
#include <vector>
#include <map>

/**
 * Lightweight transformer/LLM for DNA sequence processing
 * Simplified attention mechanism for sequence analysis
 */
class LightweightLLM {
public:
    LightweightLLM(int vocab_size = 4,  // A, T, G, C
                   int d_model = 64,
                   int num_heads = 4,
                   int num_layers = 2);
    
    /**
     * Process sequence through transformer
     * @param sequence Input DNA sequence
     * @return Encoded representation
     */
    std::vector<double> encode(const std::string& sequence);
    
    /**
     * Multi-head self-attention
     * @param input Input embeddings
     * @return Attention output
     */
    std::vector<std::vector<double>> selfAttention(const std::vector<std::vector<double>>& input);
    
    /**
     * Search for pattern using transformer embeddings
     * @param sequence DNA sequence
     * @param pattern Pattern to find
     * @param threshold Similarity threshold
     * @return Positions and similarity scores
     */
    std::vector<std::pair<size_t, double>> searchPattern(const std::string& sequence,
                                                        const std::string& pattern,
                                                        double threshold = 0.7);
    
    /**
     * Generate sequence embedding
     * @param sequence Input sequence
     * @return Fixed-size embedding vector
     */
    std::vector<double> generateEmbedding(const std::string& sequence);
    
private:
    int vocab_size_;
    int d_model_;
    int num_heads_;
    int num_layers_;
    
    // Token embeddings
    std::map<char, std::vector<double>> token_embeddings_;
    
    // Position embeddings
    std::vector<std::vector<double>> position_embeddings_;
    
    // Attention weights (simplified)
    std::vector<std::vector<std::vector<double>>> attention_weights_;
    
    /**
     * Tokenize DNA sequence
     */
    std::vector<char> tokenize(const std::string& sequence);
    
    /**
     * Embed tokens
     */
    std::vector<std::vector<double>> embedTokens(const std::vector<char>& tokens);
    
    /**
     * Add position encoding
     */
    void addPositionEncoding(std::vector<std::vector<double>>& embeddings);
    
    /**
     * Initialize embeddings
     */
    void initializeEmbeddings();
    
    /**
     * Scaled dot-product attention
     */
    std::vector<std::vector<double>> scaledDotProductAttention(
        const std::vector<std::vector<double>>& Q,
        const std::vector<std::vector<double>>& K,
        const std::vector<std::vector<double>>& V);
    
    /**
     * Layer normalization (simplified)
     */
    void layerNorm(std::vector<std::vector<double>>& input);
    
    /**
     * Feed-forward network
     */
    std::vector<std::vector<double>> feedForward(const std::vector<std::vector<double>>& input);
};

#endif // LIGHTWEIGHT_LLM_H

