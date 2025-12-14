#include "sequence_generator.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <cctype>
#include <random>

SequenceGenerator::SequenceGenerator(unsigned int seed) 
    : rng(seed == 0 ? std::random_device{}() : seed), base_dist(0, 3) {
}

std::string SequenceGenerator::generateHighEntropy(int length) {
    std::string sequence;
    sequence.reserve(length);
    
    const char bases[] = {'A', 'T', 'G', 'C'};
    
    for (int i = 0; i < length; ++i) {
        sequence += bases[base_dist(rng)];
    }
    
    return sequence;
}

std::string SequenceGenerator::generateLowEntropy(int length, int pattern_length) {
    std::string sequence;
    sequence.reserve(length);
    
    // Generate a random pattern
    std::string pattern = generateHighEntropy(pattern_length);
    
    // Repeat the pattern
    int num_full_repeats = length / pattern_length;
    for (int i = 0; i < num_full_repeats; ++i) {
        sequence += pattern;
    }
    
    // Add partial repeat if needed
    int remainder = length % pattern_length;
    if (remainder > 0) {
        sequence += pattern.substr(0, remainder);
    }
    
    return sequence;
}

std::string SequenceGenerator::generateModerateComplexity(int length, double repetition_ratio) {
    std::string sequence;
    sequence.reserve(length);
    
    int repetitive_length = static_cast<int>(length * repetition_ratio);
    int random_length = length - repetitive_length;
    
    // Add repetitive part
    if (repetitive_length > 0) {
        sequence += generateLowEntropy(repetitive_length, 4);
    }
    
    // Add random part
    if (random_length > 0) {
        sequence += generateHighEntropy(random_length);
    }
    
    return sequence;
}

std::string SequenceGenerator::generateTandemRepeat(const std::string& repeat_unit, int num_repeats) {
    std::string sequence;
    sequence.reserve(repeat_unit.length() * num_repeats);
    
    for (int i = 0; i < num_repeats; ++i) {
        sequence += repeat_unit;
    }
    
    return sequence;
}

std::string SequenceGenerator::generateWithGCContent(int length, double gc_content) {
    std::string sequence;
    sequence.reserve(length);
    
    int gc_count = static_cast<int>(length * gc_content);
    int at_count = length - gc_count;
    
    // Create pool of bases
    std::string bases;
    for (int i = 0; i < gc_count; ++i) {
        bases += randomGCBase();
    }
    for (int i = 0; i < at_count; ++i) {
        bases += randomATBase();
    }
    
    // Shuffle
    std::shuffle(bases.begin(), bases.end(), rng);
    
    return bases;
}

double SequenceGenerator::calculateEntropy(const std::string& sequence) {
    if (sequence.empty()) {
        return 0.0;
    }
    
    std::map<char, int> counts;
    for (char c : sequence) {
        char upper = std::toupper(c);
        if (upper == 'A' || upper == 'T' || upper == 'G' || upper == 'C') {
            counts[upper]++;
        }
    }
    
    double entropy = 0.0;
    double length = static_cast<double>(sequence.length());
    
    for (const auto& pair : counts) {
        double probability = static_cast<double>(pair.second) / length;
        if (probability > 0.0) {
            entropy -= probability * std::log2(probability);
        }
    }
    
    return entropy;
}

double SequenceGenerator::calculateGCContent(const std::string& sequence) {
    if (sequence.empty()) {
        return 0.0;
    }
    
    int gc_count = 0;
    int total = 0;
    
    for (char c : sequence) {
        char upper = std::toupper(c);
        if (upper == 'G' || upper == 'C') {
            gc_count++;
        }
        if (upper == 'A' || upper == 'T' || upper == 'G' || upper == 'C') {
            total++;
        }
    }
    
    if (total == 0) {
        return 0.0;
    }
    
    return static_cast<double>(gc_count) / total;
}

char SequenceGenerator::randomBase() {
    const char bases[] = {'A', 'T', 'G', 'C'};
    return bases[base_dist(rng)];
}

char SequenceGenerator::randomGCBase() {
    const char gc_bases[] = {'G', 'C'};
    std::uniform_int_distribution<int> gc_dist(0, 1);
    return gc_bases[gc_dist(rng)];
}

char SequenceGenerator::randomATBase() {
    const char at_bases[] = {'A', 'T'};
    std::uniform_int_distribution<int> at_dist(0, 1);
    return at_bases[at_dist(rng)];
}

