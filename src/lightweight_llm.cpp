#include "lightweight_llm.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <cctype>

LightweightLLM::LightweightLLM(int vocab_size, int d_model, int num_heads, int num_layers)
    : vocab_size_(vocab_size), d_model_(d_model), num_heads_(num_heads), num_layers_(num_layers) {
    initializeEmbeddings();
}

std::vector<double> LightweightLLM::encode(const std::string& sequence) {
    std::vector<char> tokens = tokenize(sequence);
    std::vector<std::vector<double>> embeddings = embedTokens(tokens);
    
    addPositionEncoding(embeddings);
    
    // Apply transformer layers
    for (int layer = 0; layer < num_layers_; ++layer) {
        // Self-attention
        std::vector<std::vector<double>> attn_output = selfAttention(embeddings);
        
        // Residual connection (simplified)
        for (size_t i = 0; i < embeddings.size() && i < attn_output.size(); ++i) {
            for (size_t j = 0; j < embeddings[i].size() && j < attn_output[i].size(); ++j) {
                embeddings[i][j] += attn_output[i][j];
            }
        }
        
        layerNorm(embeddings);
        
        // Feed-forward
        std::vector<std::vector<double>> ff_output = feedForward(embeddings);
        
        // Residual connection
        for (size_t i = 0; i < embeddings.size() && i < ff_output.size(); ++i) {
            for (size_t j = 0; j < embeddings[i].size() && j < ff_output[i].size(); ++j) {
                embeddings[i][j] += ff_output[i][j];
            }
        }
        
        layerNorm(embeddings);
    }
    
    // Pool to fixed-size embedding (mean pooling)
    std::vector<double> final_embedding(d_model_, 0.0);
    if (!embeddings.empty()) {
        for (const auto& emb : embeddings) {
            for (size_t i = 0; i < emb.size() && i < final_embedding.size(); ++i) {
                final_embedding[i] += emb[i];
            }
        }
        for (double& val : final_embedding) {
            val /= embeddings.size();
        }
    }
    
    return final_embedding;
}

std::vector<std::vector<double>> LightweightLLM::selfAttention(const std::vector<std::vector<double>>& input) {
    // Simplified self-attention: use input as Q, K, V
    return scaledDotProductAttention(input, input, input);
}

std::vector<std::pair<size_t, double>> LightweightLLM::searchPattern(const std::string& sequence,
                                                                    const std::string& pattern,
                                                                    double threshold) {
    std::vector<std::pair<size_t, double>> results;
    
    if (pattern.length() > sequence.length()) {
        return results;
    }
    
    // Generate pattern embedding
    std::vector<double> pattern_emb = generateEmbedding(pattern);
    
    // Slide window over sequence
    for (size_t i = 0; i <= sequence.length() - pattern.length(); ++i) {
        std::string window = sequence.substr(i, pattern.length());
        std::vector<double> window_emb = generateEmbedding(window);
        
        // Calculate cosine similarity
        double similarity = 0.0;
        double norm1 = 0.0, norm2 = 0.0;
        for (size_t j = 0; j < pattern_emb.size() && j < window_emb.size(); ++j) {
            similarity += pattern_emb[j] * window_emb[j];
            norm1 += pattern_emb[j] * pattern_emb[j];
            norm2 += window_emb[j] * window_emb[j];
        }
        
        double denominator = std::sqrt(norm1) * std::sqrt(norm2);
        if (denominator > 0) {
            similarity /= denominator;
        }
        
        if (similarity >= threshold) {
            results.push_back({i, similarity});
        }
    }
    
    return results;
}

std::vector<double> LightweightLLM::generateEmbedding(const std::string& sequence) {
    return encode(sequence);
}

std::vector<char> LightweightLLM::tokenize(const std::string& sequence) {
    std::vector<char> tokens;
    for (char c : sequence) {
        char base = std::toupper(c);
        if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
            tokens.push_back(base);
        }
    }
    return tokens;
}

std::vector<std::vector<double>> LightweightLLM::embedTokens(const std::vector<char>& tokens) {
    std::vector<std::vector<double>> embeddings;
    
    for (char token : tokens) {
        if (token_embeddings_.find(token) != token_embeddings_.end()) {
            embeddings.push_back(token_embeddings_[token]);
        } else {
            // Unknown token - use zero vector
            embeddings.push_back(std::vector<double>(d_model_, 0.0));
        }
    }
    
    return embeddings;
}

void LightweightLLM::addPositionEncoding(std::vector<std::vector<double>>& embeddings) {
    for (size_t pos = 0; pos < embeddings.size() && pos < position_embeddings_.size(); ++pos) {
        for (size_t i = 0; i < embeddings[pos].size() && i < position_embeddings_[pos].size(); ++i) {
            embeddings[pos][i] += position_embeddings_[pos][i];
        }
    }
}

void LightweightLLM::initializeEmbeddings() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-0.1, 0.1);
    
    // Initialize token embeddings for A, T, G, C
    char bases[] = {'A', 'T', 'G', 'C'};
    for (char base : bases) {
        std::vector<double> embedding(d_model_);
        for (int i = 0; i < d_model_; ++i) {
            embedding[i] = dist(gen);
        }
        token_embeddings_[base] = embedding;
    }
    
    // Initialize position embeddings (for up to 1000 positions)
    int max_pos = 1000;
    position_embeddings_.resize(max_pos);
    for (int pos = 0; pos < max_pos; ++pos) {
        position_embeddings_[pos].resize(d_model_);
        for (int i = 0; i < d_model_; ++i) {
            if (i % 2 == 0) {
                position_embeddings_[pos][i] = std::sin(pos / std::pow(10000.0, 2.0 * i / d_model_));
            } else {
                position_embeddings_[pos][i] = std::cos(pos / std::pow(10000.0, 2.0 * i / d_model_));
            }
        }
    }
}

std::vector<std::vector<double>> LightweightLLM::scaledDotProductAttention(
    const std::vector<std::vector<double>>& Q,
    const std::vector<std::vector<double>>& K,
    const std::vector<std::vector<double>>& V) {
    
    int seq_len = Q.size();
    int d_k = (Q.empty()) ? d_model_ : Q[0].size();
    
    std::vector<std::vector<double>> scores(seq_len, std::vector<double>(seq_len, 0.0));
    
    // Compute attention scores: Q * K^T / sqrt(d_k)
    for (int i = 0; i < seq_len; ++i) {
        for (int j = 0; j < seq_len; ++j) {
            double score = 0.0;
            for (size_t k = 0; k < Q[i].size() && k < K[j].size(); ++k) {
                score += Q[i][k] * K[j][k];
            }
            scores[i][j] = score / std::sqrt(static_cast<double>(d_k));
        }
    }
    
    // Apply softmax
    for (int i = 0; i < seq_len; ++i) {
        double max_score = *std::max_element(scores[i].begin(), scores[i].end());
        double sum = 0.0;
        for (int j = 0; j < seq_len; ++j) {
            scores[i][j] = std::exp(scores[i][j] - max_score);
            sum += scores[i][j];
        }
        for (int j = 0; j < seq_len; ++j) {
            scores[i][j] /= sum;
        }
    }
    
    // Apply to values: scores * V
    std::vector<std::vector<double>> output(seq_len, std::vector<double>(d_k, 0.0));
    for (int i = 0; i < seq_len; ++i) {
        for (int j = 0; j < seq_len; ++j) {
            for (size_t k = 0; k < V[j].size() && k < output[i].size(); ++k) {
                output[i][k] += scores[i][j] * V[j][k];
            }
        }
    }
    
    return output;
}

void LightweightLLM::layerNorm(std::vector<std::vector<double>>& input) {
    // Simplified layer normalization
    for (auto& vec : input) {
        double mean = 0.0;
        for (double val : vec) {
            mean += val;
        }
        mean /= vec.size();
        
        double variance = 0.0;
        for (double val : vec) {
            variance += (val - mean) * (val - mean);
        }
        variance /= vec.size();
        double std_dev = std::sqrt(variance + 1e-8);
        
        for (double& val : vec) {
            val = (val - mean) / std_dev;
        }
    }
}

std::vector<std::vector<double>> LightweightLLM::feedForward(const std::vector<std::vector<double>>& input) {
    // Simplified feed-forward: linear transformation
    std::vector<std::vector<double>> output = input;
    
    // In real implementation, would have weight matrices
    // Here we just return input (placeholder)
    
    return output;
}

