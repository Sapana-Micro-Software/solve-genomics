#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <fstream>
#include <sstream>
#include "ExactMatch.h"
#include "NaiveSearch.h"
#include "FuzzySearch.h"
#include "SmithWaterman.h"
#include "NeedlemanWunsch.h"
#include "sequence_generator.h"
#include "utils.h"

using namespace std;
using namespace DNAUtils;

struct BenchmarkResult {
    string algorithm_name;
    string sequence_type;
    int sequence_length;
    int pattern_length;
    double execution_time_ms;
    int64_t matches_found;
    double entropy;
    double gc_content;
    bool success;
    
    BenchmarkResult() : execution_time_ms(0.0), matches_found(0), 
                       entropy(0.0), gc_content(0.0), success(false) {}
};

class BenchmarkRunner {
private:
    vector<BenchmarkResult> results;
    SequenceGenerator generator;
    
public:
    BenchmarkRunner() : generator(42) {}  // Fixed seed for reproducibility
    
    void runSearchBenchmark(const string& algorithm_name,
                           const string& sequence,
                           const string& pattern,
                           const string& sequence_type) {
        BenchmarkResult result;
        result.algorithm_name = algorithm_name;
        result.sequence_type = sequence_type;
        result.sequence_length = sequence.length();
        result.pattern_length = pattern.length();
        result.entropy = SequenceGenerator::calculateEntropy(sequence);
        result.gc_content = SequenceGenerator::calculateGCContent(sequence);
        
        auto start = chrono::high_resolution_clock::now();
        
        if (algorithm_name == "ExactMatch") {
            ExactMatch matcher;
            SearchResult search_result = matcher.search(sequence, pattern);
            result.matches_found = search_result.count;
            result.success = true;
        } else if (algorithm_name == "NaiveSearch") {
            NaiveSearch searcher;
            SearchResult search_result = searcher.search(sequence, pattern);
            result.matches_found = search_result.count;
            result.success = true;
        } else if (algorithm_name == "FuzzySearch_0") {
            FuzzySearch fuzzy;
            FuzzySearchResult fuzzy_result = fuzzy.search(sequence, pattern, 0);
            result.matches_found = fuzzy_result.count;
            result.success = true;
        } else if (algorithm_name == "FuzzySearch_1") {
            FuzzySearch fuzzy;
            FuzzySearchResult fuzzy_result = fuzzy.search(sequence, pattern, 1);
            result.matches_found = fuzzy_result.count;
            result.success = true;
        } else if (algorithm_name == "FuzzySearch_2") {
            FuzzySearch fuzzy;
            FuzzySearchResult fuzzy_result = fuzzy.search(sequence, pattern, 2);
            result.matches_found = fuzzy_result.count;
            result.success = true;
        }
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        result.execution_time_ms = duration.count() / 1000.0;
        
        results.push_back(result);
    }
    
    void runAlignmentBenchmark(const string& algorithm_name,
                              const string& seq1,
                              const string& seq2,
                              const string& sequence_type) {
        BenchmarkResult result;
        result.algorithm_name = algorithm_name;
        result.sequence_type = sequence_type;
        result.sequence_length = seq1.length();
        result.pattern_length = seq2.length();
        result.entropy = SequenceGenerator::calculateEntropy(seq1);
        result.gc_content = SequenceGenerator::calculateGCContent(seq1);
        
        auto start = chrono::high_resolution_clock::now();
        
        if (algorithm_name == "SmithWaterman") {
            SmithWaterman aligner;
            aligner.setScoring(2, -1, -1);
            AlignmentResult align_result = aligner.align(seq1, seq2);
            result.matches_found = (align_result.score > 0) ? 1 : 0;
            result.success = true;
        } else if (algorithm_name == "NeedlemanWunsch") {
            NeedlemanWunsch aligner;
            aligner.setScoring(2, -1, -1);
            AlignmentResult align_result = aligner.align(seq1, seq2);
            result.matches_found = (align_result.score > 0) ? 1 : 0;
            result.success = true;
        }
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        result.execution_time_ms = duration.count() / 1000.0;
        
        results.push_back(result);
    }
    
    void printResults() {
        cout << "\n" << string(100, '=') << endl;
        cout << "BENCHMARK RESULTS" << endl;
        cout << string(100, '=') << endl;
        
        cout << left << setw(20) << "Algorithm" 
             << setw(20) << "Sequence Type"
             << setw(12) << "Seq Length"
             << setw(12) << "Pat Length"
             << setw(15) << "Time (ms)"
             << setw(12) << "Matches"
             << setw(12) << "Entropy"
             << setw(12) << "GC Content"
             << endl;
        cout << string(100, '-') << endl;
        
        for (const auto& result : results) {
            cout << left << setw(20) << result.algorithm_name
                 << setw(20) << result.sequence_type
                 << setw(12) << result.sequence_length
                 << setw(12) << result.pattern_length
                 << fixed << setprecision(3) << setw(15) << result.execution_time_ms
                 << setw(12) << result.matches_found
                 << setprecision(2) << setw(12) << result.entropy
                 << setprecision(2) << setw(12) << result.gc_content
                 << endl;
        }
    }
    
    void saveResultsToCSV(const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Could not open file " << filename << endl;
            return;
        }
        
        // Header
        file << "Algorithm,SequenceType,SequenceLength,PatternLength,Time_ms,Matches,Entropy,GCContent\n";
        
        // Data
        for (const auto& result : results) {
            file << result.algorithm_name << ","
                 << result.sequence_type << ","
                 << result.sequence_length << ","
                 << result.pattern_length << ","
                 << fixed << setprecision(6) << result.execution_time_ms << ","
                 << result.matches_found << ","
                 << setprecision(4) << result.entropy << ","
                 << setprecision(4) << result.gc_content << "\n";
        }
        
        file.close();
        cout << "\nResults saved to " << filename << endl;
    }
    
    void clearResults() {
        results.clear();
    }
};

void runComplexityBenchmarks() {
    BenchmarkRunner runner;
    SequenceGenerator gen(42);
    
    cout << "Running Complexity Benchmarks..." << endl;
    cout << "=================================" << endl;
    
    // Test different sequence lengths
    vector<int> sequence_lengths = {100, 500, 1000, 5000, 10000};
    vector<int> pattern_lengths = {4, 8, 16, 32};
    
    for (int seq_len : sequence_lengths) {
        for (int pat_len : pattern_lengths) {
            if (pat_len > seq_len) continue;
            
            // High entropy sequences
            string high_entropy_seq = gen.generateHighEntropy(seq_len);
            string high_entropy_pattern = gen.generateHighEntropy(pat_len);
            
            // Low entropy (repetitive) sequences
            string low_entropy_seq = gen.generateLowEntropy(seq_len, 4);
            string low_entropy_pattern = "ATCG";  // Simple repetitive pattern
            
            // Moderate complexity
            string moderate_seq = gen.generateModerateComplexity(seq_len, 0.5);
            string moderate_pattern = gen.generateHighEntropy(pat_len);
            
            cout << "\nTesting: Seq=" << seq_len << ", Pat=" << pat_len << endl;
            
            // Test search algorithms
            vector<string> search_algorithms = {"ExactMatch", "NaiveSearch", 
                                               "FuzzySearch_0", "FuzzySearch_1", "FuzzySearch_2"};
            
            for (const string& algo : search_algorithms) {
                runner.runSearchBenchmark(algo, high_entropy_seq, high_entropy_pattern, "HighEntropy");
                runner.runSearchBenchmark(algo, low_entropy_seq, low_entropy_pattern, "LowEntropy");
                runner.runSearchBenchmark(algo, moderate_seq, moderate_pattern, "Moderate");
            }
            
            // Test alignment algorithms (use smaller sequences for alignment)
            if (seq_len <= 1000) {
                string seq1_high = gen.generateHighEntropy(seq_len);
                string seq2_high = gen.generateHighEntropy(seq_len);
                string seq1_low = gen.generateLowEntropy(seq_len, 4);
                string seq2_low = gen.generateLowEntropy(seq_len, 4);
                
                runner.runAlignmentBenchmark("SmithWaterman", seq1_high, seq2_high, "HighEntropy");
                runner.runAlignmentBenchmark("SmithWaterman", seq1_low, seq2_low, "LowEntropy");
                runner.runAlignmentBenchmark("NeedlemanWunsch", seq1_high, seq2_high, "HighEntropy");
                runner.runAlignmentBenchmark("NeedlemanWunsch", seq1_low, seq2_low, "LowEntropy");
            }
        }
    }
    
    runner.printResults();
    runner.saveResultsToCSV("benchmark_results.csv");
}

void runDetailedBenchmark() {
    BenchmarkRunner runner;
    SequenceGenerator gen(42);
    
    cout << "\nRunning Detailed Benchmark..." << endl;
    cout << "==============================" << endl;
    
    // Fixed size for detailed analysis
    int seq_len = 1000;
    int pat_len = 10;
    
    // Generate sequences
    string high_entropy = gen.generateHighEntropy(seq_len);
    string low_entropy = gen.generateLowEntropy(seq_len, 4);
    string tandem_repeat = gen.generateTandemRepeat("ATCG", seq_len / 4);
    string pattern = gen.generateHighEntropy(pat_len);
    
    cout << "\nSequence Characteristics:" << endl;
    cout << "High Entropy - Entropy: " << fixed << setprecision(4) 
         << SequenceGenerator::calculateEntropy(high_entropy) << " bits" << endl;
    cout << "Low Entropy - Entropy: " << SequenceGenerator::calculateEntropy(low_entropy) << " bits" << endl;
    cout << "Tandem Repeat - Entropy: " << SequenceGenerator::calculateEntropy(tandem_repeat) << " bits" << endl;
    
    // Run multiple iterations for statistical significance
    int iterations = 10;
    
    cout << "\nRunning " << iterations << " iterations per algorithm..." << endl;
    
    for (int i = 0; i < iterations; ++i) {
        // High entropy tests
        runner.runSearchBenchmark("ExactMatch", high_entropy, pattern, "HighEntropy");
        runner.runSearchBenchmark("NaiveSearch", high_entropy, pattern, "HighEntropy");
        runner.runSearchBenchmark("FuzzySearch_1", high_entropy, pattern, "HighEntropy");
        
        // Low entropy tests
        runner.runSearchBenchmark("ExactMatch", low_entropy, "ATCG", "LowEntropy");
        runner.runSearchBenchmark("NaiveSearch", low_entropy, "ATCG", "LowEntropy");
        runner.runSearchBenchmark("FuzzySearch_1", low_entropy, "ATCG", "LowEntropy");
        
        // Tandem repeat tests
        runner.runSearchBenchmark("ExactMatch", tandem_repeat, "ATCG", "TandemRepeat");
        runner.runSearchBenchmark("NaiveSearch", tandem_repeat, "ATCG", "TandemRepeat");
    }
    
    runner.printResults();
}

int main(int argc, char* argv[]) {
    cout << "DNA Sequence Alignment Benchmark Suite" << endl;
    cout << "=======================================" << endl;
    
    if (argc > 1 && string(argv[1]) == "detailed") {
        runDetailedBenchmark();
    } else {
        runComplexityBenchmarks();
    }
    
    return 0;
}

