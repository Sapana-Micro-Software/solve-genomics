#ifndef EMBEDDING_SEARCH_H
#define EMBEDDING_SEARCH_H

#include <string>
#include <vector>
#include <map>

/**
 * Vector embedding for DNA sequences
 * Converts sequences to fixed-size vectors for similarity search
 */
class SequenceEmbedding {
public:
    /**
     * Generate embedding for a DNA sequence
     * @param sequence Input DNA sequence
     * @param embedding_size Size of embedding vector
     * @return Embedding vector
     */
    static std::vector<double> generateEmbedding(const std::string& sequence, int embedding_size = 64);
    
    /**
     * Generate k-mer based embedding (count-based)
     * @param sequence Input sequence
     * @param k K-mer size (default: 3)
     * @return Embedding vector (4^k dimensions)
     */
    static std::vector<double> generateKmerEmbedding(const std::string& sequence, int k = 3);
    
    /**
     * Generate frequency-based embedding
     * @param sequence Input sequence
     * @return Embedding with base frequencies and patterns
     */
    static std::vector<double> generateFrequencyEmbedding(const std::string& sequence);
    
    /**
     * Calculate cosine similarity between two embeddings
     */
    static double cosineSimilarity(const std::vector<double>& emb1, const std::vector<double>& emb2);
    
    /**
     * Calculate Euclidean distance between embeddings
     */
    static double euclideanDistance(const std::vector<double>& emb1, const std::vector<double>& emb2);
};

/**
 * Embedding-based sequence search
 */
class EmbeddingSearch {
public:
    EmbeddingSearch(int embedding_size = 64);
    
    /**
     * Add sequence to search index
     */
    void addSequence(const std::string& sequence, size_t id);
    
    /**
     * Build index from sequences
     */
    void buildIndex();
    
    /**
     * Search for similar sequences using embeddings
     * @param query Query sequence
     * @param top_k Number of results to return
     * @return Vector of (id, similarity_score) pairs
     */
    std::vector<std::pair<size_t, double>> search(const std::string& query, int top_k = 10);
    
    /**
     * Search with similarity threshold
     * @param query Query sequence
     * @param threshold Minimum similarity threshold
     * @return Vector of (id, similarity_score) pairs
     */
    std::vector<std::pair<size_t, double>> searchThreshold(const std::string& query, double threshold);
    
private:
    int embedding_size_;
    std::map<size_t, std::string> sequences_;
    std::map<size_t, std::vector<double>> embeddings_;
    
    /**
     * Normalize embedding vector
     */
    void normalizeEmbedding(std::vector<double>& embedding);
};

#endif // EMBEDDING_SEARCH_H

