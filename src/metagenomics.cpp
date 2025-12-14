#include "metagenomics.h"
#include "SmithWaterman.h"
#include "FuzzySearch.h"
#include <algorithm>
#include <cmath>
#include <sstream>

Metagenomics::Metagenomics() {
}

double Metagenomics::calculateSimilarity(const std::string& read, const std::string& reference) {
    if (read.empty() || reference.empty()) {
        return 0.0;
    }
    
    // Use Smith-Waterman for local alignment
    SmithWaterman aligner;
    aligner.setScoring(2, -1, -1);
    AlignmentResult alignment = aligner.align(read, reference);
    
    // Calculate similarity as normalized score
    int max_score = static_cast<int>(std::min(read.length(), reference.length())) * 2;
    if (max_score == 0) {
        return 0.0;
    }
    
    return std::max(0.0, std::min(1.0, static_cast<double>(alignment.score) / max_score));
}

Metagenomics::MetagenomicResult Metagenomics::classifyReads(
    const std::vector<std::string>& reads,
    const std::map<std::string, std::string>& reference_db,
    double similarity_threshold) {
    
    MetagenomicResult result;
    result.total_reads = static_cast<int>(reads.size());
    
    std::map<std::string, int> organism_counts;
    std::map<std::string, double> organism_similarities;
    std::map<size_t, std::string> assignments;
    
    // Classify each read
    for (size_t i = 0; i < reads.size(); ++i) {
        std::string best_organism = findBestMatch(reads[i], reference_db, similarity_threshold);
        
        if (!best_organism.empty()) {
            assignments[i] = best_organism;
            organism_counts[best_organism]++;
            result.classified_reads++;
            
            // Calculate average similarity
            double sim = calculateSimilarity(reads[i], reference_db.at(best_organism));
            if (organism_similarities.find(best_organism) == organism_similarities.end()) {
                organism_similarities[best_organism] = 0.0;
            }
            organism_similarities[best_organism] += sim;
        }
    }
    
    // Build organism matches
    for (const auto& [organism, count] : organism_counts) {
        OrganismMatch match;
        match.organism_name = organism;
        match.read_count = count;
        match.abundance = static_cast<double>(count) / result.total_reads;
        
        if (organism_similarities.find(organism) != organism_similarities.end()) {
            match.similarity = organism_similarities[organism] / count;
        }
        
        result.organisms.push_back(match);
    }
    
    // Sort by abundance
    std::sort(result.organisms.begin(), result.organisms.end(),
        [](const OrganismMatch& a, const OrganismMatch& b) {
            return a.abundance > b.abundance;
        });
    
    // Calculate diversity
    std::map<std::string, double> abundances;
    for (const auto& org : result.organisms) {
        abundances[org.organism_name] = org.abundance;
    }
    result.diversity_index = calculateShannonDiversity(abundances);
    
    return result;
}

std::map<std::string, double> Metagenomics::estimateAbundance(
    const std::vector<std::string>& reads,
    const std::map<size_t, std::string>& organism_assignments) {
    
    std::map<std::string, int> counts;
    
    for (const auto& [read_idx, organism] : organism_assignments) {
        counts[organism]++;
    }
    
    std::map<std::string, double> abundances;
    int total = static_cast<int>(organism_assignments.size());
    
    if (total > 0) {
        for (const auto& [organism, count] : counts) {
            abundances[organism] = static_cast<double>(count) / total;
        }
    }
    
    return abundances;
}

double Metagenomics::calculateShannonDiversity(const std::map<std::string, double>& abundances) {
    double diversity = 0.0;
    
    for (const auto& [organism, abundance] : abundances) {
        if (abundance > 0.0) {
            diversity -= abundance * std::log2(abundance);
        }
    }
    
    return diversity;
}

double Metagenomics::calculateSimpsonDiversity(const std::map<std::string, double>& abundances) {
    double simpson = 0.0;
    
    for (const auto& [organism, abundance] : abundances) {
        simpson += abundance * abundance;
    }
    
    return 1.0 - simpson;  // Simpson diversity index
}

std::map<std::string, std::vector<std::string>> Metagenomics::buildTaxonomicTree(
    const std::vector<OrganismMatch>& organisms) {
    
    std::map<std::string, std::vector<std::string>> tree;
    Metagenomics metagenomics;
    
    for (const auto& org : organisms) {
        if (!org.taxonomy.empty()) {
            std::vector<std::string> levels = metagenomics.parseTaxonomy(org.taxonomy);
            
            for (size_t i = 0; i < levels.size(); ++i) {
                std::string level = levels[i];
                if (i > 0) {
                    std::string parent = levels[i-1];
                    if (tree.find(parent) == tree.end()) {
                        tree[parent] = {};
                    }
                    // Add child if not already present
                    bool found = false;
                    for (const std::string& child : tree[parent]) {
                        if (child == level) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        tree[parent].push_back(level);
                    }
                }
            }
        }
    }
    
    return tree;
}

std::vector<std::pair<std::string, size_t>> Metagenomics::findMarkerGenes(
    const std::string& sequence,
    const std::map<std::string, std::string>& marker_patterns) {
    
    std::vector<std::pair<std::string, size_t>> markers;
    
    FuzzySearch fuzzy;
    
    for (const auto& [marker_name, pattern] : marker_patterns) {
        FuzzySearch fuzzy_searcher;
        FuzzySearchResult fuzzy_result = fuzzy_searcher.search(sequence, pattern, 2);  // Allow 2 errors
        SearchResult result;
        result.positions = fuzzy_result.positions;
        result.count = fuzzy_result.count;
        
        for (int pos : result.positions) {
            markers.push_back({marker_name, static_cast<size_t>(pos)});
        }
    }
    
    return markers;
}

Metagenomics::MetagenomicResult Metagenomics::assignTaxonomy(
    const std::vector<std::string>& reads,
    const std::map<std::string, std::string>& marker_db) {
    
    MetagenomicResult result;
    result.total_reads = static_cast<int>(reads.size());
    
    std::map<std::string, int> taxonomy_counts;
    std::map<std::string, std::vector<std::string>> taxonomy_reads;
    
    for (size_t i = 0; i < reads.size(); ++i) {
        std::vector<std::pair<std::string, size_t>> markers = findMarkerGenes(reads[i], marker_db);
        
        if (!markers.empty()) {
            // Use first marker found
            std::string taxonomy = marker_db.at(markers[0].first);
            taxonomy_counts[taxonomy]++;
            taxonomy_reads[taxonomy].push_back(reads[i]);
            result.classified_reads++;
        }
    }
    
    // Build organism matches from taxonomy
    for (const auto& [taxonomy, count] : taxonomy_counts) {
        OrganismMatch match;
        match.taxonomy = taxonomy;
        match.read_count = count;
        match.abundance = static_cast<double>(count) / result.total_reads;
        match.similarity = 0.8;  // Default similarity for marker-based assignment
        
        result.organisms.push_back(match);
        result.taxonomy_counts[taxonomy] = count;
    }
    
    // Calculate diversity
    std::map<std::string, double> abundances;
    for (const auto& org : result.organisms) {
        abundances[org.taxonomy] = org.abundance;
    }
    result.diversity_index = calculateShannonDiversity(abundances);
    
    return result;
}

std::vector<std::string> Metagenomics::parseTaxonomy(const std::string& taxonomy) {
    std::vector<std::string> levels;
    std::stringstream ss(taxonomy);
    std::string level;
    
    while (std::getline(ss, level, ';')) {
        if (!level.empty()) {
            levels.push_back(level);
        }
    }
    
    return levels;
}

std::string Metagenomics::findBestMatch(const std::string& read,
                                       const std::map<std::string, std::string>& reference_db,
                                       double threshold) {
    std::string best_organism;
    double best_similarity = 0.0;
    
    for (const auto& [organism, reference] : reference_db) {
        double similarity = calculateSimilarity(read, reference);
        
        if (similarity > best_similarity && similarity >= threshold) {
            best_similarity = similarity;
            best_organism = organism;
        }
    }
    
    return best_organism;
}

