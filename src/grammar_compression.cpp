#include "grammar_compression.h"
#include <algorithm>
#include <sstream>
#include <iomanip>

GrammarCompression::GrammarCompression(int min_pattern_length, int min_occurrences)
    : min_pattern_length_(min_pattern_length), min_occurrences_(min_occurrences) {
}

CompressionResult GrammarCompression::compress(const std::string& sequence) {
    CompressionResult result;
    result.is_lossless = true;
    result.original_size = sequence.length();
    
    if (sequence.empty()) {
        result.compressed_data = "";
        result.compressed_size = 0;
        result.compression_ratio = 1.0;
        return result;
    }
    
    // Find repeating patterns
    std::map<std::string, std::vector<size_t>> patterns = 
        findRepeatingPatterns(sequence, min_pattern_length_, min_occurrences_);
    
    if (patterns.empty()) {
        // No repeating patterns found, return original
        result.compressed_data = sequence;
        result.compressed_size = sequence.length();
        result.compression_ratio = 1.0;
        return result;
    }
    
    // Build grammar rules
    result.grammar = buildGrammarRules(sequence, patterns);
    
    // Create pattern to rule mapping
    std::map<std::string, std::string> pattern_to_rule;
    int rule_num = 1;
    for (const auto& pattern : patterns) {
        pattern_to_rule[pattern.first] = generateRuleName(rule_num++);
    }
    
    // Replace patterns with rule references
    result.compressed_data = replacePatternsWithRules(sequence, pattern_to_rule);
    
    // Store dictionary (pattern -> rule mapping)
    for (const auto& entry : pattern_to_rule) {
        result.dictionary[entry.second] = entry.first;
    }
    
    // Calculate compressed size (compressed data + grammar overhead)
    size_t grammar_size = 0;
    for (const auto& rule : result.grammar) {
        grammar_size += rule.length();
    }
    result.compressed_size = result.compressed_data.length() + grammar_size;
    result.compression_ratio = calculateCompressionRatio(
        result.original_size, result.compressed_size);
    
    return result;
}

std::string GrammarCompression::decompress(const CompressionResult& compressed) {
    if (compressed.compressed_data.empty()) {
        return "";
    }
    
    std::string decompressed = compressed.compressed_data;
    
    // Replace rule references with actual patterns
    // Process rules in reverse order to handle nested rules
    for (auto it = compressed.grammar.rbegin(); it != compressed.grammar.rend(); ++it) {
        const std::string& rule = *it;
        
        // Parse rule: "R1 -> ATCG"
        size_t arrow_pos = rule.find(" -> ");
        if (arrow_pos != std::string::npos) {
            std::string rule_name = rule.substr(0, arrow_pos);
            std::string rule_value = rule.substr(arrow_pos + 4);
            
            // Replace all occurrences of rule_name with rule_value
            size_t pos = 0;
            while ((pos = decompressed.find(rule_name, pos)) != std::string::npos) {
                decompressed.replace(pos, rule_name.length(), rule_value);
                pos += rule_value.length();
            }
        }
    }
    
    return decompressed;
}

std::map<std::string, std::vector<size_t>> GrammarCompression::findRepeatingPatterns(
    const std::string& sequence,
    int min_length,
    int min_count) {
    
    std::map<std::string, std::vector<size_t>> patterns;
    
    // Find all substrings of length min_length and longer
    for (size_t len = min_length; len <= sequence.length() / min_count; ++len) {
        std::map<std::string, std::vector<size_t>> current_patterns;
        
        // Slide window of length len
        for (size_t i = 0; i <= sequence.length() - len; ++i) {
            std::string pattern = sequence.substr(i, len);
            current_patterns[pattern].push_back(i);
        }
        
        // Keep patterns that occur at least min_count times
        for (const auto& entry : current_patterns) {
            if (entry.second.size() >= static_cast<size_t>(min_count)) {
                // Check if this pattern is not a substring of an already found pattern
                bool is_substring = false;
                for (const auto& existing : patterns) {
                    if (existing.first.find(entry.first) != std::string::npos &&
                        existing.first != entry.first) {
                        is_substring = true;
                        break;
                    }
                }
                
                if (!is_substring) {
                    patterns[entry.first] = entry.second;
                }
            }
        }
    }
    
    return patterns;
}

std::vector<std::string> GrammarCompression::buildGrammarRules(
    const std::string& /* sequence */,
    const std::map<std::string, std::vector<size_t>>& patterns) {
    
    std::vector<std::string> grammar;
    int rule_num = 1;
    
    // Sort patterns by length (longer first) and frequency
    std::vector<std::pair<std::string, size_t>> sorted_patterns;
    for (const auto& pattern : patterns) {
        sorted_patterns.push_back({pattern.first, pattern.second.size()});
    }
    
    std::sort(sorted_patterns.begin(), sorted_patterns.end(),
        [](const std::pair<std::string, size_t>& a, 
           const std::pair<std::string, size_t>& b) {
            if (a.second != b.second) {
                return a.second > b.second;  // More frequent first
            }
            return a.first.length() > b.first.length();  // Longer first
        });
    
    // Create grammar rules
    for (const auto& pattern_entry : sorted_patterns) {
        std::string rule_name = generateRuleName(rule_num++);
        std::string rule = rule_name + " -> " + pattern_entry.first;
        grammar.push_back(rule);
    }
    
    return grammar;
}

std::string GrammarCompression::replacePatternsWithRules(
    const std::string& sequence,
    const std::map<std::string, std::string>& pattern_to_rule) {
    
    std::string result = sequence;
    
    // Sort patterns by length (longer first) to avoid partial replacements
    std::vector<std::pair<std::string, std::string>> sorted_rules;
    for (const auto& entry : pattern_to_rule) {
        sorted_rules.push_back({entry.first, entry.second});
    }
    
    std::sort(sorted_rules.begin(), sorted_rules.end(),
        [](const std::pair<std::string, std::string>& a,
           const std::pair<std::string, std::string>& b) {
            return a.first.length() > b.first.length();
        });
    
    // Replace patterns (longer first)
    for (const auto& rule : sorted_rules) {
        const std::string& pattern = rule.first;
        const std::string& rule_name = rule.second;
        
        size_t pos = 0;
        while ((pos = result.find(pattern, pos)) != std::string::npos) {
            result.replace(pos, pattern.length(), rule_name);
            pos += rule_name.length();
        }
    }
    
    return result;
}

std::string GrammarCompression::generateRuleName(int rule_number) {
    std::ostringstream oss;
    oss << "R" << rule_number;
    return oss.str();
}

