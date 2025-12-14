#ifndef SEQUENCE_GENERATOR_H
#define SEQUENCE_GENERATOR_H

#include <string>
#include <random>

/**
 * Generator for DNA sequences with different complexity characteristics
 */
class SequenceGenerator {
public:
    SequenceGenerator(unsigned int seed = 0);
    
    /**
     * Generate high entropy (random) DNA sequence
     * High entropy = more random, less repetitive, more information content
     * @param length Length of sequence to generate
     * @return Random DNA sequence
     */
    std::string generateHighEntropy(int length);
    
    /**
     * Generate low entropy (highly repetitive) DNA sequence
     * Low entropy = repetitive patterns, less information content
     * @param length Length of sequence to generate
     * @param pattern_length Length of repeating pattern
     * @return Repetitive DNA sequence
     */
    std::string generateLowEntropy(int length, int pattern_length = 4);
    
    /**
     * Generate sequence with moderate complexity
     * Mix of random and repetitive regions
     * @param length Length of sequence to generate
     * @param repetition_ratio Ratio of repetitive to random (0.0 to 1.0)
     * @return Mixed complexity sequence
     */
    std::string generateModerateComplexity(int length, double repetition_ratio = 0.5);
    
    /**
     * Generate sequence with tandem repeats
     * @param repeat_unit The unit that repeats
     * @param num_repeats Number of times to repeat
     * @return Tandem repeat sequence
     */
    std::string generateTandemRepeat(const std::string& repeat_unit, int num_repeats);
    
    /**
     * Generate sequence with specific GC content
     * @param length Length of sequence
     * @param gc_content Desired GC content (0.0 to 1.0)
     * @return Sequence with specified GC content
     */
    std::string generateWithGCContent(int length, double gc_content = 0.5);
    
    /**
     * Calculate entropy of a sequence
     * @param sequence DNA sequence
     * @return Shannon entropy in bits
     */
    static double calculateEntropy(const std::string& sequence);
    
    /**
     * Calculate GC content of a sequence
     * @param sequence DNA sequence
     * @return GC content ratio (0.0 to 1.0)
     */
    static double calculateGCContent(const std::string& sequence);
    
private:
    std::mt19937 rng;
    std::uniform_int_distribution<int> base_dist;
    
    char randomBase();
    char randomGCBase();
    char randomATBase();
};

#endif // SEQUENCE_GENERATOR_H

