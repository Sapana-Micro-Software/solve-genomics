#ifndef VARIANT_CALLING_H
#define VARIANT_CALLING_H

#include "SequenceAligner.h"
#include <string>
#include <vector>
#include <map>

/**
 * Variant Calling: Identify genetic variants (SNPs, indels, etc.)
 */
class VariantCalling {
public:
    /**
     * Types of variants
     */
    enum VariantType {
        SNP,           // Single Nucleotide Polymorphism
        INSERTION,     // Insertion of one or more bases
        DELETION,      // Deletion of one or more bases
        SUBSTITUTION,  // Multi-base substitution
        COMPLEX        // Complex variant (multiple changes)
    };
    
    /**
     * Variant information
     */
    struct Variant {
        size_t position;           // Position in reference
        VariantType type;          // Type of variant
        std::string ref_allele;    // Reference allele
        std::string alt_allele;    // Alternative allele
        double quality_score;      // Quality score (0-100)
        int depth;                 // Read depth at position
        double allele_frequency;   // Allele frequency
        
        Variant() : position(0), type(SNP), quality_score(0.0), 
                   depth(0), allele_frequency(0.0) {}
    };
    
    /**
     * Variant calling result
     */
    struct VariantResult {
        std::vector<Variant> variants;
        int total_variants;
        int snp_count;
        int indel_count;
        double average_quality;
        
        VariantResult() : total_variants(0), snp_count(0), 
                         indel_count(0), average_quality(0.0) {}
    };
    
    VariantCalling();
    
    /**
     * Call variants by comparing sequence to reference
     * @param reference Reference sequence
     * @param sequence Query sequence
     * @param min_quality Minimum quality score threshold
     * @return VariantResult with all detected variants
     */
    VariantResult callVariants(const std::string& reference,
                               const std::string& sequence,
                               double min_quality = 20.0);
    
    /**
     * Call variants from multiple reads (pileup)
     * @param reference Reference sequence
     * @param reads Vector of read sequences
     * @param min_quality Minimum quality score
     * @param min_depth Minimum read depth
     * @return VariantResult with consensus variants
     */
    VariantResult callVariantsPileup(const std::string& reference,
                                    const std::vector<std::string>& reads,
                                    double min_quality = 20.0,
                                    int min_depth = 3);
    
    /**
     * Classify variant type
     * @param ref_allele Reference allele
     * @param alt_allele Alternative allele
     * @return VariantType
     */
    static VariantType classifyVariant(const std::string& ref_allele,
                                      const std::string& alt_allele);
    
    /**
     * Calculate quality score for variant
     * @param reference Reference sequence
     * @param sequence Query sequence
     * @param position Variant position
     * @param variant Variant information
     * @return Quality score (0-100)
     */
    static double calculateQualityScore(const std::string& reference,
                                       const std::string& sequence,
                                       size_t position,
                                       const Variant& variant);
    
    /**
     * Filter variants by quality and depth
     * @param variants Input variants
     * @param min_quality Minimum quality
     * @param min_depth Minimum depth
     * @return Filtered variants
     */
    static std::vector<Variant> filterVariants(const std::vector<Variant>& variants,
                                               double min_quality,
                                               int min_depth);
    
    /**
     * Annotate variant (predict effect)
     * @param variant Variant to annotate
     * @param reference Reference sequence
     * @return Annotation string
     */
    static std::string annotateVariant(const Variant& variant,
                                      const std::string& reference);
    
private:
    /**
     * Align sequences and find differences
     */
    std::vector<Variant> findDifferences(const std::string& reference,
                                         const std::string& sequence);
    
    /**
     * Calculate read depth at position
     */
    int calculateDepth(const std::vector<std::string>& reads, size_t position);
    
    /**
     * Calculate allele frequency
     */
    double calculateAlleleFrequency(const std::vector<std::string>& reads,
                                    size_t position,
                                    const std::string& alt_allele);
};

#endif // VARIANT_CALLING_H

