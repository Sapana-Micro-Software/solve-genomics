#include "variant_calling.h"
#include "NeedlemanWunsch.h"
#include <algorithm>
#include <cmath>
#include <cctype>

VariantCalling::VariantCalling() {
}

VariantCalling::VariantResult VariantCalling::callVariants(const std::string& reference,
                                                           const std::string& sequence,
                                                           double min_quality) {
    VariantResult result;
    
    if (reference.empty() || sequence.empty()) {
        return result;
    }
    
    // Align sequences
    NeedlemanWunsch aligner;
    aligner.setScoring(2, -1, -1);
    AlignmentResult alignment = aligner.align(reference, sequence);
    
    // Find differences
    std::vector<Variant> variants = findDifferences(reference, sequence);
    
    // Calculate quality scores and filter
    for (auto& variant : variants) {
        variant.quality_score = calculateQualityScore(reference, sequence, 
                                                     variant.position, variant);
        variant.depth = 1;  // Single read
        variant.allele_frequency = 1.0;
    }
    
    // Filter by quality
    variants = filterVariants(variants, min_quality, 1);
    
    result.variants = variants;
    result.total_variants = static_cast<int>(variants.size());
    
    for (const auto& variant : variants) {
        if (variant.type == SNP) {
            result.snp_count++;
        } else if (variant.type == INSERTION || variant.type == DELETION) {
            result.indel_count++;
        }
        result.average_quality += variant.quality_score;
    }
    
    if (result.total_variants > 0) {
        result.average_quality /= result.total_variants;
    }
    
    return result;
}

VariantCalling::VariantResult VariantCalling::callVariantsPileup(
    const std::string& reference,
    const std::vector<std::string>& reads,
    double min_quality,
    int min_depth) {
    
    VariantResult result;
    
    if (reference.empty() || reads.empty()) {
        return result;
    }
    
    // Build consensus from reads
    std::map<size_t, std::map<char, int>> pileup;
    
    for (const std::string& read : reads) {
        // Simple alignment (in practice, would use proper alignment)
        for (size_t i = 0; i < read.length() && i < reference.length(); ++i) {
            char base = std::toupper(read[i]);
            if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
                pileup[i][base]++;
            }
        }
    }
    
    // Call variants from pileup
    std::vector<Variant> variants;
    
    for (const auto& [pos, bases] : pileup) {
        if (pos >= reference.length()) continue;
        
        char ref_base = std::toupper(reference[pos]);
        int total_depth = 0;
        char most_common = ref_base;
        int max_count = 0;
        
        for (const auto& [base, count] : bases) {
            total_depth += count;
            if (count > max_count) {
                max_count = count;
                most_common = base;
            }
        }
        
        if (total_depth < min_depth) continue;
        
        // Check if variant exists
        if (most_common != ref_base && max_count > 0) {
            Variant variant;
            variant.position = pos;
            variant.ref_allele = std::string(1, ref_base);
            variant.alt_allele = std::string(1, most_common);
            variant.type = classifyVariant(variant.ref_allele, variant.alt_allele);
            variant.depth = total_depth;
            variant.allele_frequency = static_cast<double>(max_count) / total_depth;
            variant.quality_score = calculateQualityScore(reference, reads[0], pos, variant);
            
            if (variant.quality_score >= min_quality && variant.depth >= min_depth) {
                variants.push_back(variant);
            }
        }
    }
    
    result.variants = variants;
    result.total_variants = static_cast<int>(variants.size());
    
    for (const auto& variant : variants) {
        if (variant.type == SNP) {
            result.snp_count++;
        } else if (variant.type == INSERTION || variant.type == DELETION) {
            result.indel_count++;
        }
        result.average_quality += variant.quality_score;
    }
    
    if (result.total_variants > 0) {
        result.average_quality /= result.total_variants;
    }
    
    return result;
}

VariantCalling::VariantType VariantCalling::classifyVariant(const std::string& ref_allele,
                                                            const std::string& alt_allele) {
    if (ref_allele.empty() || alt_allele.empty()) {
        return COMPLEX;
    }
    
    if (ref_allele.length() == 1 && alt_allele.length() == 1) {
        return SNP;
    } else if (ref_allele.length() < alt_allele.length()) {
        return INSERTION;
    } else if (ref_allele.length() > alt_allele.length()) {
        return DELETION;
    } else if (ref_allele.length() == alt_allele.length() && ref_allele.length() > 1) {
        return SUBSTITUTION;
    }
    
    return COMPLEX;
}

double VariantCalling::calculateQualityScore(const std::string& reference,
                                             const std::string& sequence,
                                             size_t position,
                                             const Variant& variant) {
    if (position >= reference.length() || position >= sequence.length()) {
        return 0.0;
    }
    
    // Phred-like quality score
    double score = 0.0;
    
    // Base quality (simplified - would use actual quality scores in practice)
    score += 30.0;
    
    // Depth penalty
    if (variant.depth < 5) {
        score -= (5 - variant.depth) * 2.0;
    }
    
    // Allele frequency penalty
    if (variant.allele_frequency < 0.2) {
        score -= (0.2 - variant.allele_frequency) * 20.0;
    }
    
    // Context quality (check surrounding bases)
    int context_matches = 0;
    int context_size = 3;
    for (int i = -context_size; i <= context_size; ++i) {
        if (i == 0) continue;
        size_t ref_pos = position + i;
        size_t seq_pos = position + i;
        
        if (ref_pos < reference.length() && seq_pos < sequence.length()) {
            if (std::toupper(reference[ref_pos]) == std::toupper(sequence[seq_pos])) {
                context_matches++;
            }
        }
    }
    
    score += (context_matches / (2.0 * context_size)) * 10.0;
    
    return std::max(0.0, std::min(100.0, score));
}

std::vector<VariantCalling::Variant> VariantCalling::filterVariants(
    const std::vector<Variant>& variants,
    double min_quality,
    int min_depth) {
    
    std::vector<Variant> filtered;
    
    for (const auto& variant : variants) {
        if (variant.quality_score >= min_quality && variant.depth >= min_depth) {
            filtered.push_back(variant);
        }
    }
    
    return filtered;
}

std::string VariantCalling::annotateVariant(const Variant& variant,
                                           const std::string& reference) {
    std::string annotation;
    
    switch (variant.type) {
        case SNP:
            annotation = "SNP: " + variant.ref_allele + " -> " + variant.alt_allele;
            break;
        case INSERTION:
            annotation = "INSERTION: +" + std::to_string(variant.alt_allele.length() - 
                        variant.ref_allele.length()) + "bp";
            break;
        case DELETION:
            annotation = "DELETION: -" + std::to_string(variant.ref_allele.length() - 
                        variant.alt_allele.length()) + "bp";
            break;
        case SUBSTITUTION:
            annotation = "SUBSTITUTION: " + variant.ref_allele + " -> " + variant.alt_allele;
            break;
        case COMPLEX:
            annotation = "COMPLEX: " + variant.ref_allele + " -> " + variant.alt_allele;
            break;
    }
    
    annotation += " at position " + std::to_string(variant.position);
    annotation += " (Quality: " + std::to_string(static_cast<int>(variant.quality_score)) + ")";
    
    return annotation;
}

std::vector<VariantCalling::Variant> VariantCalling::findDifferences(
    const std::string& reference,
    const std::string& sequence) {
    
    std::vector<Variant> variants;
    
    // Simple difference finding (would use proper alignment in practice)
    size_t min_len = std::min(reference.length(), sequence.length());
    
    for (size_t i = 0; i < min_len; ++i) {
        char ref_base = std::toupper(reference[i]);
        char seq_base = std::toupper(sequence[i]);
        
        if (ref_base != seq_base && 
            (ref_base == 'A' || ref_base == 'T' || ref_base == 'G' || ref_base == 'C') &&
            (seq_base == 'A' || seq_base == 'T' || seq_base == 'G' || seq_base == 'C')) {
            
            Variant variant;
            variant.position = i;
            variant.ref_allele = std::string(1, ref_base);
            variant.alt_allele = std::string(1, seq_base);
            variant.type = classifyVariant(variant.ref_allele, variant.alt_allele);
            
            variants.push_back(variant);
        }
    }
    
    // Handle indels (simplified)
    if (reference.length() != sequence.length()) {
        size_t max_len = std::max(reference.length(), sequence.length());
        if (max_len > min_len) {
            Variant variant;
            variant.position = min_len;
            
            if (reference.length() > sequence.length()) {
                variant.type = DELETION;
                variant.ref_allele = reference.substr(min_len);
                variant.alt_allele = "";
            } else {
                variant.type = INSERTION;
                variant.ref_allele = "";
                variant.alt_allele = sequence.substr(min_len);
            }
            
            variants.push_back(variant);
        }
    }
    
    return variants;
}

int VariantCalling::calculateDepth(const std::vector<std::string>& reads, size_t position) {
    int depth = 0;
    
    for (const std::string& read : reads) {
        if (position < read.length()) {
            char base = std::toupper(read[position]);
            if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
                depth++;
            }
        }
    }
    
    return depth;
}

double VariantCalling::calculateAlleleFrequency(const std::vector<std::string>& reads,
                                                size_t position,
                                                const std::string& alt_allele) {
    if (reads.empty() || alt_allele.empty()) {
        return 0.0;
    }
    
    int alt_count = 0;
    int total_count = 0;
    
    for (const std::string& read : reads) {
        if (position < read.length()) {
            char base = std::toupper(read[position]);
            if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
                total_count++;
                if (base == std::toupper(alt_allele[0])) {
                    alt_count++;
                }
            }
        }
    }
    
    if (total_count == 0) {
        return 0.0;
    }
    
    return static_cast<double>(alt_count) / total_count;
}

