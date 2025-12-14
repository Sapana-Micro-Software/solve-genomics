#ifndef VECTOR_SIMILARITY_H
#define VECTOR_SIMILARITY_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <cmath>

/**
 * Vector-based similarity and distance algorithms for DNA sequences
 */
class VectorSimilarity {
public:
    /**
     * k-Nearest Neighbors (kNN) for sequence similarity search
     */
    class KNN {
    public:
        struct Neighbor {
            size_t index;
            double distance;
            std::string sequence;
            
            Neighbor(size_t idx, double dist, const std::string& seq)
                : index(idx), distance(dist), sequence(seq) {}
            
            bool operator>(const Neighbor& other) const {
                return distance > other.distance;
            }
        };
        
        KNN(int k = 5);
        
        /**
         * Build kNN index from sequences
         * @param sequences Sequences to index
         */
        void buildIndex(const std::vector<std::string>& sequences);
        
        /**
         * Find k nearest neighbors
         * @param query Query sequence
         * @param k Number of neighbors (overrides constructor)
         * @return Vector of k nearest neighbors
         */
        std::vector<Neighbor> findNeighbors(const std::string& query, int k = -1);
        
        /**
         * Find neighbors within radius
         * @param query Query sequence
         * @param radius Maximum distance
         * @return Vector of neighbors within radius
         */
        std::vector<Neighbor> findNeighborsRadius(const std::string& query, double radius);
        
    private:
        int k_;
        std::vector<std::string> sequences_;
        std::vector<std::vector<double>> embeddings_;
        
        /**
         * Convert sequence to vector embedding
         */
        std::vector<double> sequenceToVector(const std::string& sequence);
        
        /**
         * Calculate Euclidean distance
         */
        double euclideanDistance(const std::vector<double>& v1, const std::vector<double>& v2);
    };
    
    /**
     * Approximate Nearest Neighbors (ANN) using locality-sensitive hashing
     */
    class ApproximateNearestNeighbors {
    public:
        ApproximateNearestNeighbors(int num_tables = 10, int num_bits = 8);
        
        /**
         * Build ANN index
         * @param sequences Sequences to index
         */
        void buildIndex(const std::vector<std::string>& sequences);
        
        /**
         * Find approximate nearest neighbors
         * @param query Query sequence
         * @param k Number of neighbors
         * @return Vector of (index, distance) pairs
         */
        std::vector<std::pair<size_t, double>> findNeighbors(const std::string& query, int k = 10);
        
    private:
        int num_tables_;
        int num_bits_;
        std::vector<std::string> sequences_;
        std::vector<std::vector<double>> embeddings_;
        std::vector<std::map<size_t, std::vector<size_t>>> hash_tables_;  // Table -> hash -> indices
        
        /**
         * Hash function for LSH
         */
        size_t hashLSH(const std::vector<double>& vector, int table_id);
        
        /**
         * Convert sequence to vector
         */
        std::vector<double> sequenceToVector(const std::string& sequence);
    };
    
    /**
     * Dynamic Time Warping (DTW) for sequence alignment
     */
    class DynamicTimeWarping {
    public:
        /**
         * Calculate DTW distance between two sequences
         * @param seq1 First sequence (as vector)
         * @param seq2 Second sequence (as vector)
         * @return DTW distance
         */
        static double distance(const std::vector<double>& seq1, const std::vector<double>& seq2);
        
        /**
         * Calculate DTW distance between DNA sequences
         * @param seq1 First DNA sequence
         * @param seq2 Second DNA sequence
         * @return DTW distance
         */
        static double distance(const std::string& seq1, const std::string& seq2);
        
        /**
         * Find optimal warping path
         * @param seq1 First sequence
         * @param seq2 Second sequence
         * @return Warping path (pairs of indices)
         */
        static std::vector<std::pair<int, int>> findWarpingPath(const std::string& seq1, const std::string& seq2);
        
    private:
        /**
         * Cost function for DTW
         */
        static double cost(double a, double b);
        
        /**
         * Convert DNA sequence to numeric vector
         */
        static std::vector<double> sequenceToVector(const std::string& sequence);
    };
    
    /**
     * Edit Distance in Vector Space
     */
    class VectorEditDistance {
    public:
        /**
         * Calculate edit distance between vector sequences
         * @param vec1 First vector sequence
         * @param vec2 Second vector sequence
         * @return Edit distance
         */
        static double distance(const std::vector<std::vector<double>>& vec1,
                              const std::vector<std::vector<double>>& vec2);
        
        /**
         * Calculate edit distance between DNA sequences as vectors
         * @param seq1 First DNA sequence
         * @param seq2 Second DNA sequence
         * @return Edit distance
         */
        static double distance(const std::string& seq1, const std::string& seq2);
        
        /**
         * Find alignment path
         * @param vec1 First vector sequence
         * @param vec2 Second vector sequence
         * @return Alignment path
         */
        static std::vector<std::pair<int, int>> findAlignmentPath(
            const std::vector<std::vector<double>>& vec1,
            const std::vector<std::vector<double>>& vec2);
        
    private:
        /**
         * Distance between two vectors
         */
        static double vectorDistance(const std::vector<double>& v1, const std::vector<double>& v2);
        
        /**
         * Convert DNA sequence to vector sequence
         */
        static std::vector<std::vector<double>> sequenceToVectors(const std::string& sequence);
    };
    
    /**
     * Levenshtein Distance in Vector Space
     */
    class VectorLevenshtein {
    public:
        /**
         * Calculate Levenshtein distance between vector sequences
         * @param vec1 First vector sequence
         * @param vec2 Second vector sequence
         * @param threshold Threshold for early termination
         * @return Levenshtein distance
         */
        static double distance(const std::vector<std::vector<double>>& vec1,
                              const std::vector<std::vector<double>>& vec2,
                              double threshold = -1.0);
        
        /**
         * Calculate Levenshtein distance between DNA sequences as vectors
         * @param seq1 First DNA sequence
         * @param seq2 Second DNA sequence
         * @return Levenshtein distance
         */
        static double distance(const std::string& seq1, const std::string& seq2);
        
        /**
         * Bounded Levenshtein distance
         * @param vec1 First vector sequence
         * @param vec2 Second vector sequence
         * @param max_distance Maximum distance to compute
         * @return Distance or max_distance+1 if exceeded
         */
        static double boundedDistance(const std::vector<std::vector<double>>& vec1,
                                     const std::vector<std::vector<double>>& vec2,
                                     double max_distance);
        
    private:
        /**
         * Distance between two vectors (for substitution cost)
         */
        static double vectorDistance(const std::vector<double>& v1, const std::vector<double>& v2);
        
        /**
         * Convert DNA sequence to vector sequence
         */
        static std::vector<std::vector<double>> sequenceToVectors(const std::string& sequence);
    };
};

#endif // VECTOR_SIMILARITY_H

