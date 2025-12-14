#ifndef DDMCMC_SEARCH_H
#define DDMCMC_SEARCH_H

#include "embedding_search.h"
#include <string>
#include <vector>
#include <random>
#include <map>

/**
 * Data-Driven MCMC (DDMCMC) for vector embedding search
 * Uses data distribution to guide MCMC sampling in embedding space
 */
class DDMCMCSearch {
public:
    struct DDMCMCResult {
        std::vector<size_t> matched_ids;    // IDs of matched sequences
        std::vector<double> similarities;    // Similarity scores
        std::vector<double> final_embedding; // Final embedding vector
        int iterations;                     // MCMC iterations
        double acceptance_rate;             // Proposal acceptance rate
        
        DDMCMCResult() : iterations(0), acceptance_rate(0.0) {}
    };
    
    DDMCMCSearch(int embedding_size = 64,
                int max_iterations = 1000,
                double step_size = 0.1,
                double temperature = 1.0);
    
    /**
     * Build data-driven proposal distribution from sequence embeddings
     * @param sequences Sequences to build distribution from
     */
    void buildProposalDistribution(const std::vector<std::string>& sequences);
    
    /**
     * Search using DDMCMC in embedding space
     * @param query Query sequence
     * @param top_k Number of results to return
     * @return DDMCMCResult with matched sequences
     */
    DDMCMCResult search(const std::string& query, int top_k = 10);
    
    /**
     * Propose new embedding using data-driven distribution
     * @param current_embedding Current embedding vector
     * @return Proposed embedding vector
     */
    std::vector<double> proposeEmbedding(const std::vector<double>& current_embedding);
    
    /**
     * Calculate likelihood of embedding given data distribution
     */
    double calculateLikelihood(const std::vector<double>& embedding);
    
    /**
     * Metropolis-Hastings acceptance for embeddings
     */
    bool acceptProposal(const std::vector<double>& current,
                       const std::vector<double>& proposed,
                       double current_likelihood,
                       double proposed_likelihood);
    
private:
    int embedding_size_;
    int max_iterations_;
    double step_size_;
    double temperature_;
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_dist_;
    std::normal_distribution<double> normal_dist_;
    
    // Data-driven proposal distribution
    std::vector<std::vector<double>> data_embeddings_;
    std::vector<double> mean_embedding_;
    std::vector<std::vector<double>> covariance_matrix_;  // Simplified: diagonal
    
    /**
     * Calculate mean embedding from data
     */
    void calculateMeanEmbedding();
    
    /**
     * Calculate covariance (simplified: diagonal)
     */
    void calculateCovariance();
    
    /**
     * Sample from data-driven distribution
     */
    std::vector<double> sampleFromDistribution();
    
    /**
     * Calculate cosine similarity
     */
    double cosineSimilarity(const std::vector<double>& v1, const std::vector<double>& v2);
};

#endif // DDMCMC_SEARCH_H

