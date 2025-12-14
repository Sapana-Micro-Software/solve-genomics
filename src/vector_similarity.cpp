#include "vector_similarity.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <cctype>
#include <climits>

// KNN Implementation

VectorSimilarity::KNN::KNN(int k) : k_(k) {
    if (k_ < 1) k_ = 5;
}

std::vector<double> VectorSimilarity::KNN::sequenceToVector(const std::string& sequence) {
    std::vector<double> vector(256, 0.0);  // 256-dimensional embedding
    
    // K-mer frequency encoding
    int k = 3;
    for (size_t i = 0; i <= sequence.length() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        for (char& c : kmer) {
            c = std::toupper(c);
        }
        
        // Hash k-mer to vector dimension
        size_t hash = std::hash<std::string>{}(kmer);
        vector[hash % 256] += 1.0;
    }
    
    // Normalize
    double norm = 0.0;
    for (double val : vector) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    if (norm > 0) {
        for (double& val : vector) {
            val /= norm;
        }
    }
    
    return vector;
}

double VectorSimilarity::KNN::euclideanDistance(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        return 1e10;  // Large distance for mismatched dimensions
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    
    return std::sqrt(sum);
}

void VectorSimilarity::KNN::buildIndex(const std::vector<std::string>& sequences) {
    sequences_ = sequences;
    embeddings_.clear();
    
    for (const std::string& seq : sequences) {
        embeddings_.push_back(sequenceToVector(seq));
    }
}

std::vector<VectorSimilarity::KNN::Neighbor> VectorSimilarity::KNN::findNeighbors(
    const std::string& query, int k) {
    
    if (k < 0) k = k_;
    if (sequences_.empty()) {
        return {};
    }
    
    std::vector<double> query_vector = sequenceToVector(query);
    
    // Priority queue for k nearest (max heap)
    std::priority_queue<Neighbor, std::vector<Neighbor>, std::greater<Neighbor>> pq;
    
    for (size_t i = 0; i < sequences_.size(); ++i) {
        double dist = euclideanDistance(query_vector, embeddings_[i]);
        
        if (static_cast<int>(pq.size()) < k) {
            pq.push(Neighbor(i, dist, sequences_[i]));
        } else if (dist < pq.top().distance) {
            pq.pop();
            pq.push(Neighbor(i, dist, sequences_[i]));
        }
    }
    
    // Extract results
    std::vector<Neighbor> neighbors;
    while (!pq.empty()) {
        neighbors.push_back(pq.top());
        pq.pop();
    }
    
    std::reverse(neighbors.begin(), neighbors.end());
    return neighbors;
}

std::vector<VectorSimilarity::KNN::Neighbor> VectorSimilarity::KNN::findNeighborsRadius(
    const std::string& query, double radius) {
    
    std::vector<Neighbor> neighbors;
    std::vector<double> query_vector = sequenceToVector(query);
    
    for (size_t i = 0; i < sequences_.size(); ++i) {
        double dist = euclideanDistance(query_vector, embeddings_[i]);
        if (dist <= radius) {
            neighbors.push_back(Neighbor(i, dist, sequences_[i]));
        }
    }
    
    std::sort(neighbors.begin(), neighbors.end(),
        [](const Neighbor& a, const Neighbor& b) {
            return a.distance < b.distance;
        });
    
    return neighbors;
}

// Approximate Nearest Neighbors Implementation

VectorSimilarity::ApproximateNearestNeighbors::ApproximateNearestNeighbors(int num_tables, int num_bits)
    : num_tables_(num_tables), num_bits_(num_bits) {
    if (num_tables_ < 1) num_tables_ = 10;
    if (num_bits_ < 1) num_bits_ = 8;
}

std::vector<double> VectorSimilarity::ApproximateNearestNeighbors::sequenceToVector(const std::string& sequence) {
    std::vector<double> vector(64, 0.0);
    
    // Frequency-based encoding
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
        for (int i = 0; i < 4 && i < 64; ++i) {
            vector[i] = static_cast<double>(counts[i]) / total;
        }
    }
    
    // K-mer features
    int k = 2;
    for (size_t i = 0; i <= sequence.length() - k && i + 4 < 64; ++i) {
        std::string kmer = sequence.substr(i, k);
        size_t hash = std::hash<std::string>{}(kmer);
        vector[4 + (hash % 60)] += 1.0;
    }
    
    return vector;
}

size_t VectorSimilarity::ApproximateNearestNeighbors::hashLSH(const std::vector<double>& vector, int table_id) {
    std::mt19937 rng(static_cast<unsigned>(table_id));
    std::normal_distribution<double> dist(0.0, 1.0);
    
    size_t hash = 0;
    for (int bit = 0; bit < num_bits_; ++bit) {
        double projection = 0.0;
        for (size_t i = 0; i < vector.size(); ++i) {
            double weight = dist(rng);
            projection += vector[i] * weight;
        }
        
        hash |= (projection > 0 ? 1ULL : 0ULL) << bit;
    }
    
    return hash;
}

void VectorSimilarity::ApproximateNearestNeighbors::buildIndex(const std::vector<std::string>& sequences) {
    sequences_ = sequences;
    embeddings_.clear();
    hash_tables_.clear();
    hash_tables_.resize(num_tables_);
    
    for (const std::string& seq : sequences) {
        embeddings_.push_back(sequenceToVector(seq));
    }
    
    // Build hash tables
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (int table = 0; table < num_tables_; ++table) {
            size_t hash = hashLSH(embeddings_[i], table);
            hash_tables_[table][hash].push_back(i);
        }
    }
}

std::vector<std::pair<size_t, double>> VectorSimilarity::ApproximateNearestNeighbors::findNeighbors(
    const std::string& query, int k) {
    
    if (sequences_.empty()) {
        return {};
    }
    
    std::vector<double> query_vector = sequenceToVector(query);
    std::map<size_t, int> candidate_counts;
    
    // Collect candidates from hash tables
    for (int table = 0; table < num_tables_; ++table) {
        size_t hash = hashLSH(query_vector, table);
        if (hash_tables_[table].find(hash) != hash_tables_[table].end()) {
            for (size_t idx : hash_tables_[table][hash]) {
                candidate_counts[idx]++;
            }
        }
    }
    
    // Calculate distances for candidates
    std::vector<std::pair<size_t, double>> candidates;
    for (const auto& [idx, count] : candidate_counts) {
        double dist = 0.0;
        for (size_t i = 0; i < query_vector.size() && i < embeddings_[idx].size(); ++i) {
            double diff = query_vector[i] - embeddings_[idx][i];
            dist += diff * diff;
        }
        dist = std::sqrt(dist);
        candidates.push_back({idx, dist});
    }
    
    // Sort and return top k
    std::sort(candidates.begin(), candidates.end(),
        [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
            return a.second < b.second;
        });
    
    if (static_cast<int>(candidates.size()) > k) {
        candidates.resize(k);
    }
    
    return candidates;
}

// Dynamic Time Warping Implementation

double VectorSimilarity::DynamicTimeWarping::cost(double a, double b) {
    return std::abs(a - b);
}

std::vector<double> VectorSimilarity::DynamicTimeWarping::sequenceToVector(const std::string& sequence) {
    std::vector<double> vector;
    vector.reserve(sequence.length());
    
    for (char c : sequence) {
        char upper_c = std::toupper(c);
        if (upper_c == 'A') vector.push_back(0.0);
        else if (upper_c == 'T') vector.push_back(1.0);
        else if (upper_c == 'G') vector.push_back(2.0);
        else if (upper_c == 'C') vector.push_back(3.0);
        else vector.push_back(0.0);
    }
    
    return vector;
}

double VectorSimilarity::DynamicTimeWarping::distance(const std::vector<double>& seq1,
                                                      const std::vector<double>& seq2) {
    int m = static_cast<int>(seq1.size());
    int n = static_cast<int>(seq2.size());
    
    if (m == 0) return n;
    if (n == 0) return m;
    
    // DTW matrix
    std::vector<std::vector<double>> dtw(m + 1, std::vector<double>(n + 1, 1e10));
    dtw[0][0] = 0.0;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            double c = cost(seq1[i-1], seq2[j-1]);
            dtw[i][j] = c + std::min({dtw[i-1][j],      // Insertion
                                     dtw[i][j-1],      // Deletion
                                     dtw[i-1][j-1]});  // Match
        }
    }
    
    return dtw[m][n];
}

double VectorSimilarity::DynamicTimeWarping::distance(const std::string& seq1, const std::string& seq2) {
    std::vector<double> vec1 = sequenceToVector(seq1);
    std::vector<double> vec2 = sequenceToVector(seq2);
    return distance(vec1, vec2);
}

std::vector<std::pair<int, int>> VectorSimilarity::DynamicTimeWarping::findWarpingPath(
    const std::string& seq1, const std::string& seq2) {
    
    std::vector<double> vec1 = sequenceToVector(seq1);
    std::vector<double> vec2 = sequenceToVector(seq2);
    
    int m = static_cast<int>(vec1.size());
    int n = static_cast<int>(vec2.size());
    
    // Build DTW matrix and traceback
    std::vector<std::vector<double>> dtw(m + 1, std::vector<double>(n + 1, 1e10));
    std::vector<std::vector<int>> path_i(m + 1, std::vector<int>(n + 1, -1));
    std::vector<std::vector<int>> path_j(m + 1, std::vector<int>(n + 1, -1));
    
    dtw[0][0] = 0.0;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            double c = cost(vec1[i-1], vec2[j-1]);
            
            double min_val = dtw[i-1][j];
            int min_i = i-1, min_j = j;
            
            if (dtw[i][j-1] < min_val) {
                min_val = dtw[i][j-1];
                min_i = i; min_j = j-1;
            }
            if (dtw[i-1][j-1] < min_val) {
                min_val = dtw[i-1][j-1];
                min_i = i-1; min_j = j-1;
            }
            
            dtw[i][j] = c + min_val;
            path_i[i][j] = min_i;
            path_j[i][j] = min_j;
        }
    }
    
    // Traceback
    std::vector<std::pair<int, int>> path;
    int i = m, j = n;
    
    while (i >= 0 && j >= 0) {
        path.push_back({i-1, j-1});
        if (i == 0 && j == 0) break;
        int next_i = path_i[i][j];
        int next_j = path_j[i][j];
        i = next_i;
        j = next_j;
    }
    
    std::reverse(path.begin(), path.end());
    return path;
}

// Vector Edit Distance Implementation

double VectorSimilarity::VectorEditDistance::vectorDistance(const std::vector<double>& v1,
                                                            const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        return 1e10;
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    
    return std::sqrt(sum);
}

std::vector<std::vector<double>> VectorSimilarity::VectorEditDistance::sequenceToVectors(
    const std::string& sequence) {
    
    std::vector<std::vector<double>> vectors;
    vectors.reserve(sequence.length());
    
    for (char c : sequence) {
        std::vector<double> vec(4, 0.0);  // One-hot encoding
        char upper_c = std::toupper(c);
        if (upper_c == 'A') vec[0] = 1.0;
        else if (upper_c == 'T') vec[1] = 1.0;
        else if (upper_c == 'G') vec[2] = 1.0;
        else if (upper_c == 'C') vec[3] = 1.0;
        vectors.push_back(vec);
    }
    
    return vectors;
}

double VectorSimilarity::VectorEditDistance::distance(
    const std::vector<std::vector<double>>& vec1,
    const std::vector<std::vector<double>>& vec2) {
    
    int m = static_cast<int>(vec1.size());
    int n = static_cast<int>(vec2.size());
    
    if (m == 0) return n;
    if (n == 0) return m;
    
    std::vector<std::vector<double>> dp(m + 1, std::vector<double>(n + 1, 1e10));
    
    dp[0][0] = 0.0;
    for (int i = 1; i <= m; ++i) dp[i][0] = i;
    for (int j = 1; j <= n; ++j) dp[0][j] = j;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            double sub_cost = vectorDistance(vec1[i-1], vec2[j-1]);
            double del_cost = 1.0;
            double ins_cost = 1.0;
            
            dp[i][j] = std::min({dp[i-1][j-1] + sub_cost,
                               dp[i-1][j] + del_cost,
                               dp[i][j-1] + ins_cost});
        }
    }
    
    return dp[m][n];
}

double VectorSimilarity::VectorEditDistance::distance(const std::string& seq1, const std::string& seq2) {
    std::vector<std::vector<double>> vec1 = sequenceToVectors(seq1);
    std::vector<std::vector<double>> vec2 = sequenceToVectors(seq2);
    return distance(vec1, vec2);
}

std::vector<std::pair<int, int>> VectorSimilarity::VectorEditDistance::findAlignmentPath(
    const std::vector<std::vector<double>>& vec1,
    const std::vector<std::vector<double>>& vec2) {
    
    int m = static_cast<int>(vec1.size());
    int n = static_cast<int>(vec2.size());
    
    std::vector<std::vector<double>> dp(m + 1, std::vector<double>(n + 1, 1e10));
    std::vector<std::vector<int>> path_i(m + 1, std::vector<int>(n + 1, -1));
    std::vector<std::vector<int>> path_j(m + 1, std::vector<int>(n + 1, -1));
    
    dp[0][0] = 0.0;
    for (int i = 1; i <= m; ++i) {
        dp[i][0] = i;
        path_i[i][0] = i-1;
        path_j[i][0] = 0;
    }
    for (int j = 1; j <= n; ++j) {
        dp[0][j] = j;
        path_i[0][j] = 0;
        path_j[0][j] = j-1;
    }
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            double sub_cost = vectorDistance(vec1[i-1], vec2[j-1]);
            double del_cost = 1.0;
            double ins_cost = 1.0;
            
            double min_val = dp[i-1][j-1] + sub_cost;
            int min_i = i-1, min_j = j-1;
            
            if (dp[i-1][j] + del_cost < min_val) {
                min_val = dp[i-1][j] + del_cost;
                min_i = i-1; min_j = j;
            }
            if (dp[i][j-1] + ins_cost < min_val) {
                min_val = dp[i][j-1] + ins_cost;
                min_i = i; min_j = j-1;
            }
            
            dp[i][j] = min_val;
            path_i[i][j] = min_i;
            path_j[i][j] = min_j;
        }
    }
    
    // Traceback
    std::vector<std::pair<int, int>> path;
    int i = m, j = n;
    
    while (i >= 0 && j >= 0) {
        path.push_back({i-1, j-1});
        if (i == 0 && j == 0) break;
        int next_i = path_i[i][j];
        int next_j = path_j[i][j];
        i = next_i;
        j = next_j;
    }
    
    std::reverse(path.begin(), path.end());
    return path;
}

// Vector Levenshtein Distance Implementation

double VectorSimilarity::VectorLevenshtein::vectorDistance(const std::vector<double>& v1,
                                                           const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        return 1.0;  // Mismatch
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    
    return std::sqrt(sum) > 0.1 ? 1.0 : 0.0;  // Threshold for match
}

std::vector<std::vector<double>> VectorSimilarity::VectorLevenshtein::sequenceToVectors(
    const std::string& sequence) {
    
    std::vector<std::vector<double>> vectors;
    vectors.reserve(sequence.length());
    
    for (char c : sequence) {
        std::vector<double> vec(4, 0.0);
        char upper_c = std::toupper(c);
        if (upper_c == 'A') vec[0] = 1.0;
        else if (upper_c == 'T') vec[1] = 1.0;
        else if (upper_c == 'G') vec[2] = 1.0;
        else if (upper_c == 'C') vec[3] = 1.0;
        vectors.push_back(vec);
    }
    
    return vectors;
}

double VectorSimilarity::VectorLevenshtein::distance(
    const std::vector<std::vector<double>>& vec1,
    const std::vector<std::vector<double>>& vec2,
    double threshold) {
    
    int m = static_cast<int>(vec1.size());
    int n = static_cast<int>(vec2.size());
    
    if (m == 0) return n;
    if (n == 0) return m;
    
    // Space-optimized version
    std::vector<double> prev_row(n + 1);
    std::vector<double> curr_row(n + 1);
    
    for (int j = 0; j <= n; ++j) {
        prev_row[j] = j;
    }
    
    for (int i = 1; i <= m; ++i) {
        curr_row[0] = i;
        
        for (int j = 1; j <= n; ++j) {
            double sub_cost = vectorDistance(vec1[i-1], vec2[j-1]);
            
            curr_row[j] = std::min({prev_row[j-1] + sub_cost,  // Substitution
                                   prev_row[j] + 1.0,          // Deletion
                                   curr_row[j-1] + 1.0});      // Insertion
            
            // Early termination if threshold exceeded
            if (threshold > 0 && curr_row[j] > threshold) {
                return threshold + 1;
            }
        }
        
        prev_row = curr_row;
    }
    
    return curr_row[n];
}

double VectorSimilarity::VectorLevenshtein::distance(const std::string& seq1, const std::string& seq2) {
    std::vector<std::vector<double>> vec1 = sequenceToVectors(seq1);
    std::vector<std::vector<double>> vec2 = sequenceToVectors(seq2);
    return distance(vec1, vec2);
}

double VectorSimilarity::VectorLevenshtein::boundedDistance(
    const std::vector<std::vector<double>>& vec1,
    const std::vector<std::vector<double>>& vec2,
    double max_distance) {
    
    return distance(vec1, vec2, max_distance);
}

