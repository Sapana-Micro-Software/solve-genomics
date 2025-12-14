#include "ddmcmc_search.h"
#include <algorithm>
#include <cmath>
#include <numeric>

DDMCMCSearch::DDMCMCSearch(int embedding_size, int max_iterations, 
                          double step_size, double temperature)
    : embedding_size_(embedding_size),
      max_iterations_(max_iterations),
      step_size_(step_size),
      temperature_(temperature),
      rng_(std::random_device{}()),
      uniform_dist_(0.0, 1.0),
      normal_dist_(0.0, 1.0) {
    mean_embedding_.resize(embedding_size_, 0.0);
    covariance_matrix_.resize(embedding_size_, std::vector<double>(embedding_size_, 0.0));
}

void DDMCMCSearch::buildProposalDistribution(const std::vector<std::string>& sequences) {
    data_embeddings_.clear();
    
    // Generate embeddings for all sequences
    for (const std::string& seq : sequences) {
        std::vector<double> embedding = SequenceEmbedding::generateEmbedding(seq, embedding_size_);
        data_embeddings_.push_back(embedding);
    }
    
    if (data_embeddings_.empty()) {
        return;
    }
    
    // Calculate statistics for proposal distribution
    calculateMeanEmbedding();
    calculateCovariance();
}

DDMCMCSearch::DDMCMCResult DDMCMCSearch::search(const std::string& query, int top_k) {
    DDMCMCResult result;
    
    if (data_embeddings_.empty()) {
        return result;
    }
    
    // Start with query embedding
    std::vector<double> current_embedding = SequenceEmbedding::generateEmbedding(query, embedding_size_);
    double current_likelihood = calculateLikelihood(current_embedding);
    
    std::vector<double> best_embedding = current_embedding;
    double best_likelihood = current_likelihood;
    
    int accepted = 0;
    
    // MCMC sampling
    for (int iter = 0; iter < max_iterations_; ++iter) {
        // Propose new embedding
        std::vector<double> proposed_embedding = proposeEmbedding(current_embedding);
        double proposed_likelihood = calculateLikelihood(proposed_embedding);
        
        // Metropolis-Hastings acceptance
        if (acceptProposal(current_embedding, proposed_embedding, 
                          current_likelihood, proposed_likelihood)) {
            current_embedding = proposed_embedding;
            current_likelihood = proposed_likelihood;
            accepted++;
            
            if (current_likelihood > best_likelihood) {
                best_embedding = current_embedding;
                best_likelihood = current_likelihood;
            }
        }
    }
    
    // Find similar sequences using best embedding
    std::vector<std::pair<size_t, double>> similarities;
    for (size_t i = 0; i < data_embeddings_.size(); ++i) {
        double sim = cosineSimilarity(best_embedding, data_embeddings_[i]);
        similarities.push_back({i, sim});
    }
    
    // Sort by similarity
    std::sort(similarities.begin(), similarities.end(),
        [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
            return a.second > b.second;
        });
    
    // Return top k
    int k = std::min(top_k, static_cast<int>(similarities.size()));
    for (int i = 0; i < k; ++i) {
        result.matched_ids.push_back(similarities[i].first);
        result.similarities.push_back(similarities[i].second);
    }
    
    result.final_embedding = best_embedding;
    result.iterations = max_iterations_;
    result.acceptance_rate = static_cast<double>(accepted) / max_iterations_;
    
    return result;
}

std::vector<double> DDMCMCSearch::proposeEmbedding(const std::vector<double>& current_embedding) {
    std::vector<double> proposed = current_embedding;
    
    // Mix of random walk and data-driven proposal
    double mix_prob = uniform_dist_(rng_);
    
    if (mix_prob < 0.5 && !data_embeddings_.empty()) {
        // Data-driven proposal: sample from distribution
        std::vector<double> data_sample = sampleFromDistribution();
        for (size_t i = 0; i < proposed.size(); ++i) {
            proposed[i] = 0.7 * current_embedding[i] + 0.3 * data_sample[i];
        }
    } else {
        // Random walk proposal
        for (size_t i = 0; i < proposed.size(); ++i) {
            proposed[i] += step_size_ * normal_dist_(rng_);
        }
    }
    
    // Normalize
    double norm = 0.0;
    for (double val : proposed) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    if (norm > 0) {
        for (double& val : proposed) {
            val /= norm;
        }
    }
    
    return proposed;
}

double DDMCMCSearch::calculateLikelihood(const std::vector<double>& embedding) {
    if (data_embeddings_.empty()) {
        return 0.0;
    }
    
    // Calculate average similarity to data embeddings
    double total_similarity = 0.0;
    for (const auto& data_emb : data_embeddings_) {
        total_similarity += cosineSimilarity(embedding, data_emb);
    }
    
    return total_similarity / data_embeddings_.size();
}

bool DDMCMCSearch::acceptProposal(const std::vector<double>& current,
                                  const std::vector<double>& proposed,
                                  double current_likelihood,
                                  double proposed_likelihood) {
    if (proposed_likelihood > current_likelihood) {
        return true;
    }
    
    double delta = proposed_likelihood - current_likelihood;
    double probability = std::exp(delta / temperature_);
    
    return uniform_dist_(rng_) < probability;
}

void DDMCMCSearch::calculateMeanEmbedding() {
    if (data_embeddings_.empty()) {
        return;
    }
    
    std::fill(mean_embedding_.begin(), mean_embedding_.end(), 0.0);
    
    for (const auto& emb : data_embeddings_) {
        for (size_t i = 0; i < emb.size() && i < mean_embedding_.size(); ++i) {
            mean_embedding_[i] += emb[i];
        }
    }
    
    for (double& val : mean_embedding_) {
        val /= data_embeddings_.size();
    }
}

void DDMCMCSearch::calculateCovariance() {
    if (data_embeddings_.empty()) {
        return;
    }
    
    // Simplified: diagonal covariance
    for (size_t i = 0; i < covariance_matrix_.size(); ++i) {
        double variance = 0.0;
        for (const auto& emb : data_embeddings_) {
            if (i < emb.size()) {
                double diff = emb[i] - mean_embedding_[i];
                variance += diff * diff;
            }
        }
        variance /= data_embeddings_.size();
        covariance_matrix_[i][i] = variance + 0.01;  // Add small value for stability
    }
}

std::vector<double> DDMCMCSearch::sampleFromDistribution() {
    std::vector<double> sample(embedding_size_);
    
    for (size_t i = 0; i < sample.size(); ++i) {
        double std_dev = std::sqrt(covariance_matrix_[i][i]);
        sample[i] = mean_embedding_[i] + std_dev * normal_dist_(rng_);
    }
    
    // Normalize
    double norm = 0.0;
    for (double val : sample) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    if (norm > 0) {
        for (double& val : sample) {
            val /= norm;
        }
    }
    
    return sample;
}

double DDMCMCSearch::cosineSimilarity(const std::vector<double>& v1, const std::vector<double>& v2) {
    return SequenceEmbedding::cosineSimilarity(v1, v2);
}

