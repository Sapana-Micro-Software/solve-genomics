#ifndef METAGENOMICS_H
#define METAGENOMICS_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>

/**
 * Metagenomics: Analysis of genetic material from multiple organisms
 */
class Metagenomics {
public:
    /**
     * Organism classification result
     */
    struct OrganismMatch {
        std::string organism_name;
        std::string taxonomy;
        double similarity;
        double abundance;
        int read_count;
        
        OrganismMatch() : similarity(0.0), abundance(0.0), read_count(0) {}
    };
    
    /**
     * Metagenomic analysis result
     */
    struct MetagenomicResult {
        std::vector<OrganismMatch> organisms;
        int total_reads;
        int classified_reads;
        double diversity_index;  // Shannon diversity
        std::map<std::string, int> taxonomy_counts;
        
        MetagenomicResult() : total_reads(0), classified_reads(0), diversity_index(0.0) {}
    };
    
    Metagenomics();
    
    /**
     * Classify reads to organisms using reference database
     * @param reads Sequence reads from metagenomic sample
     * @param reference_db Map of organism names to reference sequences
     * @param similarity_threshold Minimum similarity for classification
     * @return MetagenomicResult with organism assignments
     */
    MetagenomicResult classifyReads(const std::vector<std::string>& reads,
                                   const std::map<std::string, std::string>& reference_db,
                                   double similarity_threshold = 0.8);
    
    /**
     * Estimate organism abundance
     * @param reads Reads assigned to organisms
     * @param organism_assignments Map of read index to organism name
     * @return Map of organism name to abundance (0-1)
     */
    std::map<std::string, double> estimateAbundance(
        const std::vector<std::string>& reads,
        const std::map<size_t, std::string>& organism_assignments);
    
    /**
     * Calculate diversity metrics
     * @param abundances Map of organism to abundance
     * @return Shannon diversity index
     */
    static double calculateShannonDiversity(const std::map<std::string, double>& abundances);
    
    /**
     * Calculate Simpson diversity index
     * @param abundances Map of organism to abundance
     * @return Simpson diversity index
     */
    static double calculateSimpsonDiversity(const std::map<std::string, double>& abundances);
    
    /**
     * Build taxonomic tree from classifications
     * @param organisms Organism matches
     * @return Taxonomic tree (simplified representation)
     */
    static std::map<std::string, std::vector<std::string>> buildTaxonomicTree(
        const std::vector<OrganismMatch>& organisms);
    
    /**
     * Find marker genes (16S rRNA, etc.)
     * @param sequence Sequence to search
     * @param marker_patterns Marker gene patterns
     * @return Vector of found markers
     */
    std::vector<std::pair<std::string, size_t>> findMarkerGenes(
        const std::string& sequence,
        const std::map<std::string, std::string>& marker_patterns);
    
    /**
     * Taxonomic assignment using marker genes
     * @param reads Reads to classify
     * @param marker_db Map of marker genes to taxonomy
     * @return MetagenomicResult with taxonomic assignments
     */
    MetagenomicResult assignTaxonomy(const std::vector<std::string>& reads,
                                    const std::map<std::string, std::string>& marker_db);
    
private:
    /**
     * Calculate similarity between read and reference
     */
    double calculateSimilarity(const std::string& read, const std::string& reference);
    
    /**
     * Parse taxonomy string (e.g., "Bacteria;Firmicutes;Bacillus")
     */
    std::vector<std::string> parseTaxonomy(const std::string& taxonomy);
    
    /**
     * Find best matching organism
     */
    std::string findBestMatch(const std::string& read,
                             const std::map<std::string, std::string>& reference_db,
                             double threshold);
};

#endif // METAGENOMICS_H

