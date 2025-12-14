#include "embedding_search.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <cctype>
#include <limits>

std::vector<double> SequenceEmbedding::generateEmbedding(const std::string& sequence, int embedding_size) {
    std::vector<double> embedding(embedding_size, 0.0);
    
    if (sequence.empty()) {
        return embedding;
    }
    
    // Simple hash-based embedding
    for (size_t i = 0; i < sequence.length(); ++i) {
        char base = std::toupper(sequence[i]);
        int base_val = 0;
        if (base == 'A') base_val = 0;
        else if (base == 'T') base_val = 1;
        else if (base == 'G') base_val = 2;
        else if (base == 'C') base_val = 3;
        
        // Distribute across embedding dimensions
        for (int j = 0; j < embedding_size; ++j) {
            double weight = std::sin((i * embedding_size + j) * 0.1) * (base_val + 1);
            embedding[j] += weight / sequence.length();
        }
    }
    
    // Normalize
    double norm = 0.0;
    for (double val : embedding) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    if (norm > 0) {
        for (double& val : embedding) {
            val /= norm;
        }
    }
    
    return embedding;
}

std::vector<double> SequenceEmbedding::generateKmerEmbedding(const std::string& sequence, int k) {
    int num_kmers = 1;
    for (int i = 0; i < k; ++i) {
        num_kmers *= 4;  // 4^k possible k-mers
    }
    
    std::vector<double> embedding(num_kmers, 0.0);
    
    if (sequence.length() < static_cast<size_t>(k)) {
        return embedding;
    }
    
    // Count k-mers
    for (size_t i = 0; i <= sequence.length() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        
        // Convert k-mer to index
        int index = 0;
        for (int j = 0; j < k; ++j) {
            char base = std::toupper(kmer[j]);
            int base_val = 0;
            if (base == 'A') base_val = 0;
            else if (base == 'T') base_val = 1;
            else if (base == 'G') base_val = 2;
            else if (base == 'C') base_val = 3;
            
            index = index * 4 + base_val;
        }
        
        if (index >= 0 && index < num_kmers) {
            embedding[index] += 1.0;
        }
    }
    
    // Normalize by sequence length
    double total = sequence.length() - k + 1;
    if (total > 0) {
        for (double& val : embedding) {
            val /= total;
        }
    }
    
    return embedding;
}

std::vector<double> SequenceEmbedding::generateFrequencyEmbedding(const std::string& sequence) {
    std::vector<double> embedding;
    
    if (sequence.empty()) {
        return embedding;
    }
    
    // Base frequencies
    int counts[4] = {0, 0, 0, 0};  // A, T, G, C
    for (char c : sequence) {
        char base = std::toupper(c);
        if (base == 'A') counts[0]++;
        else if (base == 'T') counts[1]++;
        else if (base == 'G') counts[2]++;
        else if (base == 'C') counts[3]++;
    }
    
    double total = sequence.length();
    for (int i = 0; i < 4; ++i) {
        embedding.push_back(counts[i] / total);
    }
    
    // GC content
    embedding.push_back((counts[2] + counts[3]) / total);
    
    // Dinucleotide frequencies
    std::map<std::string, int> dinuc_counts;
    for (size_t i = 0; i < sequence.length() - 1; ++i) {
        std::string dinuc = sequence.substr(i, 2);
        std::transform(dinuc.begin(), dinuc.end(), dinuc.begin(), ::toupper);
        dinuc_counts[dinuc]++;
    }
    
    // Add top dinucleotide frequencies
    std::vector<std::pair<std::string, int>> sorted_dinucs;
    for (const auto& entry : dinuc_counts) {
        sorted_dinucs.push_back(entry);
    }
    std::sort(sorted_dinucs.begin(), sorted_dinucs.end(),
        [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
            return a.second > b.second;
        });
    
    for (size_t i = 0; i < std::min(static_cast<size_t>(10), sorted_dinucs.size()); ++i) {
        embedding.push_back(sorted_dinucs[i].second / static_cast<double>(sequence.length() - 1));
    }
    
    return embedding;
}

double SequenceEmbedding::cosineSimilarity(const std::vector<double>& emb1, const std::vector<double>& emb2) {
    if (emb1.size() != emb2.size() || emb1.empty()) {
        return 0.0;
    }
    
    double dot_product = 0.0;
    double norm1 = 0.0;
    double norm2 = 0.0;
    
    for (size_t i = 0; i < emb1.size(); ++i) {
        dot_product += emb1[i] * emb2[i];
        norm1 += emb1[i] * emb1[i];
        norm2 += emb2[i] * emb2[i];
    }
    
    double denominator = std::sqrt(norm1) * std::sqrt(norm2);
    if (denominator == 0.0) {
        return 0.0;
    }
    
    return dot_product / denominator;
}

double SequenceEmbedding::euclideanDistance(const std::vector<double>& emb1, const std::vector<double>& emb2) {
    if (emb1.size() != emb2.size() || emb1.empty()) {
        return std::numeric_limits<double>::max();
    }
    
    double distance = 0.0;
    for (size_t i = 0; i < emb1.size(); ++i) {
        double diff = emb1[i] - emb2[i];
        distance += diff * diff;
    }
    
    return std::sqrt(distance);
}

// EmbeddingSearch implementation

EmbeddingSearch::EmbeddingSearch(int embedding_size) : embedding_size_(embedding_size) {
}

void EmbeddingSearch::addSequence(const std::string& sequence, size_t id) {
    sequences_[id] = sequence;
}

void EmbeddingSearch::buildIndex() {
    embeddings_.clear();
    
    for (const auto& entry : sequences_) {
        std::vector<double> embedding = SequenceEmbedding::generateEmbedding(entry.second, embedding_size_);
        normalizeEmbedding(embedding);
        embeddings_[entry.first] = embedding;
    }
}

std::vector<std::pair<size_t, double>> EmbeddingSearch::search(const std::string& query, int top_k) {
    std::vector<double> query_embedding = SequenceEmbedding::generateEmbedding(query, embedding_size_);
    normalizeEmbedding(query_embedding);
    
    std::vector<std::pair<size_t, double>> results;
    
    for (const auto& entry : embeddings_) {
        double similarity = SequenceEmbedding::cosineSimilarity(query_embedding, entry.second);
        results.push_back({entry.first, similarity});
    }
    
    // Sort by similarity (descending)
    std::sort(results.begin(), results.end(),
        [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
            return a.second > b.second;
        });
    
    // Return top k
    if (top_k > 0 && static_cast<size_t>(top_k) < results.size()) {
        results.resize(top_k);
    }
    
    return results;
}

std::vector<std::pair<size_t, double>> EmbeddingSearch::searchThreshold(const std::string& query, double threshold) {
    std::vector<double> query_embedding = SequenceEmbedding::generateEmbedding(query, embedding_size_);
    normalizeEmbedding(query_embedding);
    
    std::vector<std::pair<size_t, double>> results;
    
    for (const auto& entry : embeddings_) {
        double similarity = SequenceEmbedding::cosineSimilarity(query_embedding, entry.second);
        if (similarity >= threshold) {
            results.push_back({entry.first, similarity});
        }
    }
    
    // Sort by similarity
    std::sort(results.begin(), results.end(),
        [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
            return a.second > b.second;
        });
    
    return results;
}

void EmbeddingSearch::normalizeEmbedding(std::vector<double>& embedding) {
    double norm = 0.0;
    for (double val : embedding) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    
    if (norm > 0) {
        for (double& val : embedding) {
            val /= norm;
        }
    }
}

