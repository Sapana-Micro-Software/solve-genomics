#include "utils.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iomanip>

namespace DNAUtils {

    bool isValidDNA(const std::string& sequence) {
        if (sequence.empty()) {                                                               // Check if sequence is empty
            return false;                                                                     // Empty sequence is invalid
        }
        
        for (char c : sequence) {                                                             // Iterate through each character
            char upper = std::toupper(c);                                                     // Convert to uppercase for comparison
            if (upper != 'A' && upper != 'T' && upper != 'G' && upper != 'C' && upper != '-') { // Check if valid DNA base or gap
                // Allow '-' for gaps in aligned sequences
                return false;                                                                 // Invalid character found
            }
        }
        return true;                                                                          // All characters are valid
    }
    
    std::string toUpperCase(const std::string& sequence) {
        std::string result = sequence;                                                         // Create copy of sequence
        std::transform(result.begin(), result.end(), result.begin(), ::toupper);              // Transform to uppercase
        return result;                                                                        // Return uppercase sequence
    }
    
    std::string removeWhitespace(const std::string& sequence) {
        std::string result;                                                                   // Initialize result string
        result.reserve(sequence.size());                                                      // Reserve space for efficiency
        
        for (char c : sequence) {                                                             // Iterate through each character
            if (!std::isspace(c)) {                                                           // Check if not whitespace
                result += c;                                                                  // Append non-whitespace character
            }
        }
        return result;                                                                        // Return sequence without whitespace
    }
    
    std::string cleanSequence(const std::string& sequence) {
        std::string cleaned = removeWhitespace(sequence);                                     // Remove all whitespace characters
        cleaned = toUpperCase(cleaned);                                                        // Convert to uppercase
        
        // Remove gaps for validation
        std::string for_validation = cleaned;                                                  // Create copy for validation
        for_validation.erase(std::remove(for_validation.begin(), for_validation.end(), '-'), for_validation.end()); // Remove gaps
        
        if (!isValidDNA(for_validation)) {                                                    // Validate DNA sequence
            return "";                                                                        // Return empty string if invalid
        }
        
        return cleaned;                                                                       // Return cleaned sequence
    }
    
    std::string readSequenceFromFile(const std::string& filename) {
        std::ifstream file(filename);                                                         // Open file for reading
        if (!file.is_open()) {                                                                // Check if file opened successfully
            return "";                                                                        // Return empty string if file not found
        }
        
        std::stringstream buffer;                                                             // Create string stream buffer
        buffer << file.rdbuf();                                                               // Read entire file into buffer
        std::string content = buffer.str();                                                   // Convert buffer to string
        
        file.close();                                                                         // Close file
        return cleanSequence(content);                                                        // Return cleaned sequence from file
    }
    
    std::string formatAlignment(const std::string& seq1, const std::string& seq2, int line_length) {
        if (seq1.length() != seq2.length()) {
            return "Error: Aligned sequences must have the same length\n";
        }
        
        std::ostringstream result;
        int length = seq1.length();
        
        for (int i = 0; i < length; i += line_length) {
            int end = std::min(i + line_length, length);
            
            // Sequence 1
            result << "Seq1: " << seq1.substr(i, end - i) << "\n";
            
            // Match/mismatch indicator
            std::string match_line = "      ";
            for (int j = i; j < end; ++j) {
                if (seq1[j] == '-' || seq2[j] == '-') {
                    match_line += ' ';
                } else if (seq1[j] == seq2[j]) {
                    match_line += '|';
                } else {
                    match_line += '*';
                }
            }
            result << match_line << "\n";
            
            // Sequence 2
            result << "Seq2: " << seq2.substr(i, end - i) << "\n\n";
        }
        
        return result.str();
    }
    
    double calculateIdentity(const std::string& aligned_seq1, const std::string& aligned_seq2) {
        if (aligned_seq1.length() != aligned_seq2.length()) {
            return 0.0;
        }
        
        int matches = 0;
        int total = 0;
        
        for (size_t i = 0; i < aligned_seq1.length(); ++i) {
            // Skip positions with gaps in both sequences
            if (aligned_seq1[i] == '-' && aligned_seq2[i] == '-') {
                continue;
            }
            
            total++;
            if (aligned_seq1[i] == aligned_seq2[i] && aligned_seq1[i] != '-') {
                matches++;
            }
        }
        
        if (total == 0) {
            return 0.0;
        }
        
        return (static_cast<double>(matches) / total) * 100.0;
    }
}

