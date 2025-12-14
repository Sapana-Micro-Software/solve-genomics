#ifndef GRAMMAR_COMPRESSION_H
#define GRAMMAR_COMPRESSION_H

#include "compression.h"
#include <map>
#include <string>
#include <vector>

/**
 * Grammar-based compression using context-free grammar rules
 * Lossless compression that finds repeating patterns and replaces them with rules
 */
class GrammarCompression : public CompressionAlgorithm {
public:
    GrammarCompression(int min_pattern_length = 2, int min_occurrences = 2);
    
    CompressionResult compress(const std::string& sequence) override;
    std::string decompress(const CompressionResult& compressed) override;
    std::string getName() const override { return "GrammarCompression"; }
    bool isLossless() const override { return true; }
    
    /**
     * Find repeating patterns in sequence
     * @param sequence Input sequence
     * @param min_length Minimum pattern length
     * @param min_count Minimum occurrences
     * @return Map of patterns to their occurrences
     */
    std::map<std::string, std::vector<size_t>> findRepeatingPatterns(
        const std::string& sequence, 
        int min_length, 
        int min_count);
    
    /**
     * Build grammar rules from patterns
     * @param sequence Original sequence
     * @param patterns Map of patterns and their positions
     * @return Grammar rules
     */
    std::vector<std::string> buildGrammarRules(
        const std::string& sequence,
        const std::map<std::string, std::vector<size_t>>& patterns);
    
private:
    int min_pattern_length_;
    int min_occurrences_;
    
    /**
     * Replace patterns in sequence with rule references
     */
    std::string replacePatternsWithRules(
        const std::string& sequence,
        const std::map<std::string, std::string>& pattern_to_rule);
    
    /**
     * Generate rule name (e.g., R1, R2, ...)
     */
    std::string generateRuleName(int rule_number);
};

#endif // GRAMMAR_COMPRESSION_H

