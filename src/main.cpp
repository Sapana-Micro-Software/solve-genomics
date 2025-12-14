#include <iostream>
#include <iomanip>
#include <chrono>
#include "ExactMatch.h"
#include "NaiveSearch.h"
#include "SmithWaterman.h"
#include "NeedlemanWunsch.h"
#include "FuzzySearch.h"
#include "grammar_compression.h"
#include "association_list.h"
#include "lossy_compression.h"
#include "edit_distance_algorithms.h"
#include "embedding_search.h"
#include "cnn_sequence.h"
#include "lightweight_llm.h"
#include "pim_search.h"
#include "classic_string_matching.h"
#include "mcmc_search.h"
#include "ddmcmc_search.h"
#include "parallel_search.h"
#include "warp_ctc.h"
#include "aho_corasick.h"
#include "suffix_tree.h"
#include "suffix_array.h"
#include "advanced_dp.h"
#include "deep_learning_advanced.h"
#include "wu_manber.h"
#include "bitap.h"
#include "graph_based.h"
#include "heuristic_search.h"
#include "vector_similarity.h"
#include "variant_calling.h"
#include "metagenomics.h"
#include "chord_hashing.h"
#include "probabilistic_ml.h"
#include "concurrent_multi_search.h"
#include "skip_graph.h"
#include "dancing_links.h"
#include "utils.h"

using namespace std;
using namespace DNAUtils;

void printSeparator() {
    cout << string(80, '=') << endl;
}

void demonstrateExactMatch() {
    printSeparator();
    cout << "EXACT MATCH SEARCH DEMONSTRATION" << endl;
    printSeparator();
    
    string sequence = "ATGCGATCGATCGATCGATCG";
    string pattern = "ATCG";
    
    cout << "Sequence: " << sequence << endl;
    cout << "Pattern:  " << pattern << endl << endl;
    
    ExactMatch matcher;
    SearchResult result = matcher.search(sequence, pattern);
    
    cout << "Found " << result.count << " occurrence(s) at positions: ";
    for (size_t i = 0; i < result.positions.size(); ++i) {
        cout << result.positions[i];
        if (i < result.positions.size() - 1) cout << ", ";
    }
    cout << endl << endl;
}

void demonstrateNaiveSearch() {
    printSeparator();
    cout << "NAIVE SEARCH DEMONSTRATION" << endl;
    printSeparator();
    
    string sequence = "ATGCGATCGATCGATCGATCG";
    string pattern = "ATCG";
    
    cout << "Sequence: " << sequence << endl;
    cout << "Pattern:  " << pattern << endl << endl;
    
    NaiveSearch searcher;
    SearchResult result = searcher.search(sequence, pattern);
    
    cout << "Found " << result.count << " occurrence(s) at positions: ";
    for (size_t i = 0; i < result.positions.size(); ++i) {
        cout << result.positions[i];
        if (i < result.positions.size() - 1) cout << ", ";
    }
    cout << endl << endl;
}

void demonstrateSmithWaterman() {
    printSeparator();
    cout << "SMITH-WATERMAN (LOCAL ALIGNMENT) DEMONSTRATION" << endl;
    printSeparator();
    
    string seq1 = "ACGTACGTACGT";
    string seq2 = "TACGTACGTACG";
    
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl << endl;
    
    SmithWaterman aligner;
    aligner.setScoring(2, -1, -1);  // match=2, mismatch=-1, gap=-1
    
    auto start = chrono::high_resolution_clock::now();
    AlignmentResult result = aligner.align(seq1, seq2);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    cout << "Alignment Score: " << result.score << endl;
    cout << "Alignment Time: " << duration.count() << " microseconds" << endl << endl;
    
    if (!result.aligned_seq1.empty()) {
        cout << "Aligned Sequences:" << endl;
        cout << formatAlignment(result.aligned_seq1, result.aligned_seq2);
        
        double identity = calculateIdentity(result.aligned_seq1, result.aligned_seq2);
        cout << "Sequence Identity: " << fixed << setprecision(2) << identity << "%" << endl;
        cout << "Local Alignment Region:" << endl;
        cout << "  Sequence 1: positions " << result.start_pos1 << " to " << result.end_pos1 << endl;
        cout << "  Sequence 2: positions " << result.start_pos2 << " to " << result.end_pos2 << endl;
    } else {
        cout << "No significant local alignment found." << endl;
    }
    cout << endl;
}

void demonstrateNeedlemanWunsch() {
    printSeparator();
    cout << "NEEDLEMAN-WUNSCH (GLOBAL ALIGNMENT) DEMONSTRATION" << endl;
    printSeparator();
    
    string seq1 = "ACGTACGT";
    string seq2 = "ACGTACGTACGT";
    
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl << endl;
    
    NeedlemanWunsch aligner;
    aligner.setScoring(2, -1, -1);  // match=2, mismatch=-1, gap=-1
    
    auto start = chrono::high_resolution_clock::now();
    AlignmentResult result = aligner.align(seq1, seq2);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    cout << "Alignment Score: " << result.score << endl;
    cout << "Alignment Time: " << duration.count() << " microseconds" << endl << endl;
    
    if (!result.aligned_seq1.empty()) {
        cout << "Aligned Sequences:" << endl;
        cout << formatAlignment(result.aligned_seq1, result.aligned_seq2);
        
        double identity = calculateIdentity(result.aligned_seq1, result.aligned_seq2);
        cout << "Sequence Identity: " << fixed << setprecision(2) << identity << "%" << endl;
        cout << "Global Alignment:" << endl;
        cout << "  Sequence 1: full length (" << result.start_pos1 << " to " << result.end_pos1 << ")" << endl;
        cout << "  Sequence 2: full length (" << result.start_pos2 << " to " << result.end_pos2 << ")" << endl;
    }
    cout << endl;
}

void demonstrateFuzzySearch() {
    printSeparator();
    cout << "FUZZY SEARCH DEMONSTRATION" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCG";
    string pattern = "ATCC";  // One mismatch from "ATCG"
    
    cout << "Sequence: " << sequence << endl;
    cout << "Pattern:  " << pattern << endl << endl;
    
    FuzzySearch fuzzy;
    
    // Search with different distance thresholds
    for (int max_dist = 0; max_dist <= 2; ++max_dist) {
        cout << "Maximum edit distance: " << max_dist << endl;
        FuzzySearchResult result = fuzzy.search(sequence, pattern, max_dist);
        
        cout << "  Found " << result.count << " match(es)" << endl;
        if (result.count > 0) {
            cout << "  Positions and distances: ";
            for (size_t i = 0; i < result.positions.size(); ++i) {
                cout << "[" << result.positions[i] << ", dist=" << result.distances[i] << "]";
                if (i < result.positions.size() - 1) cout << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }
    
    // Demonstrate Hamming distance search
    cout << "Hamming Distance Search (same-length only):" << endl;
    string seq2 = "ATCGATCG";
    string pat2 = "ATCC";  // One substitution
    FuzzySearchResult hamming_result = fuzzy.searchHamming(seq2, pat2, 1);
    cout << "  Sequence: " << seq2 << endl;
    cout << "  Pattern:  " << pat2 << endl;
    cout << "  Found " << hamming_result.count << " match(es) with distance <= 1" << endl;
    if (hamming_result.count > 0) {
        for (size_t i = 0; i < hamming_result.positions.size(); ++i) {
            cout << "    Position " << hamming_result.positions[i] 
                 << ": distance = " << hamming_result.distances[i] << endl;
        }
    }
    cout << endl;
}

void demonstrateCompression() {
    printSeparator();
    cout << "COMPRESSION DEMONSTRATION" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCGATCGATCG";
    cout << "Original Sequence: " << sequence << endl;
    cout << "Original Size: " << sequence.length() << " bytes" << endl << endl;
    
    // Grammar-based compression (lossless)
    cout << "Grammar-Based Compression (Lossless):" << endl;
    GrammarCompression grammar_comp(2, 2);  // min_pattern=2, min_occurrences=2
    CompressionResult grammar_result = grammar_comp.compress(sequence);
    
    cout << "  Compressed: " << grammar_result.compressed_data << endl;
    cout << "  Grammar Rules:" << endl;
    for (const auto& rule : grammar_result.grammar) {
        cout << "    " << rule << endl;
    }
    cout << "  Compressed Size: " << grammar_result.compressed_size << " bytes" << endl;
    cout << "  Compression Ratio: " << fixed << setprecision(2) 
         << grammar_result.compression_ratio << endl;
    
    string decompressed = grammar_comp.decompress(grammar_result);
    cout << "  Decompressed: " << decompressed << endl;
    cout << "  Lossless: " << (decompressed == sequence ? "Yes" : "No") << endl << endl;
    
    // Lossy compression - frequency based
    cout << "Lossy Compression - Frequency Based:" << endl;
    LossyCompression lossy_comp(LossyCompression::FREQUENCY_BASED, 0.5);
    CompressionResult lossy_result = lossy_comp.compress(sequence);
    
    cout << "  Compressed: " << lossy_result.compressed_data << endl;
    cout << "  Dictionary:" << endl;
    for (const auto& entry : lossy_result.dictionary) {
        cout << "    " << entry.first << " -> " << entry.second << endl;
    }
    cout << "  Compressed Size: " << lossy_result.compressed_size << " bytes" << endl;
    cout << "  Compression Ratio: " << lossy_result.compression_ratio << endl;
    
    string lossy_decompressed = lossy_comp.decompress(lossy_result);
    cout << "  Decompressed (approximate): " << lossy_decompressed << endl;
    cout << "  Lossless: " << (lossy_decompressed == sequence ? "Yes" : "No") << endl << endl;
    
    // Association list and pattern matching
    cout << "Pattern Matching on Compressed Sequence:" << endl;
    AssociationList assoc_list;
    assoc_list.buildFromGrammar(grammar_result.grammar);
    assoc_list.buildFromDictionary(grammar_result.dictionary);
    
    CompressedPatternMatcher matcher(assoc_list);
    string pattern = "ATCG";
    vector<size_t> positions = matcher.search(grammar_result.compressed_data, pattern);
    
    cout << "  Pattern: " << pattern << endl;
    cout << "  Found at positions: ";
    for (size_t i = 0; i < positions.size(); ++i) {
        cout << positions[i];
        if (i < positions.size() - 1) cout << ", ";
    }
    cout << endl << endl;
}

void demonstrateComparison() {
    printSeparator();
    cout << "ALGORITHM COMPARISON: Local vs Global Alignment" << endl;
    printSeparator();
    
    string seq1 = "ACGTACGTACGT";
    string seq2 = "TACGTACGTACG";
    
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl << endl;
    
    // Smith-Waterman (Local)
    SmithWaterman sw;
    sw.setScoring(2, -1, -1);
    auto start_sw = chrono::high_resolution_clock::now();
    AlignmentResult result_sw = sw.align(seq1, seq2);
    auto end_sw = chrono::high_resolution_clock::now();
    auto duration_sw = chrono::duration_cast<chrono::microseconds>(end_sw - start_sw);
    
    // Needleman-Wunsch (Global)
    NeedlemanWunsch nw;
    nw.setScoring(2, -1, -1);
    auto start_nw = chrono::high_resolution_clock::now();
    AlignmentResult result_nw = nw.align(seq1, seq2);
    auto end_nw = chrono::high_resolution_clock::now();
    auto duration_nw = chrono::duration_cast<chrono::microseconds>(end_nw - start_nw);
    
    cout << "SMITH-WATERMAN (Local):" << endl;
    cout << "  Score: " << result_sw.score << endl;
    cout << "  Time: " << duration_sw.count() << " microseconds" << endl;
    if (!result_sw.aligned_seq1.empty()) {
        cout << "  Alignment length: " << result_sw.aligned_seq1.length() << " bases" << endl;
    }
    cout << endl;
    
    cout << "NEEDLEMAN-WUNSCH (Global):" << endl;
    cout << "  Score: " << result_nw.score << endl;
    cout << "  Time: " << duration_nw.count() << " microseconds" << endl;
    if (!result_nw.aligned_seq1.empty()) {
        cout << "  Alignment length: " << result_nw.aligned_seq1.length() << " bases" << endl;
    }
    cout << endl;
}

void demonstrateClassicStringMatching() {
    printSeparator();
    cout << "CLASSIC STRING MATCHING ALGORITHMS" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCGATCG";
    string pattern = "ATCG";
    
    cout << "Sequence: " << sequence << endl;
    cout << "Pattern:  " << pattern << endl << endl;
    
    // Rabin-Karp
    cout << "Rabin-Karp Algorithm:" << endl;
    auto start_rk = chrono::high_resolution_clock::now();
    SearchResult rk_result = ClassicStringMatching::RabinKarp::search(sequence, pattern);
    auto end_rk = chrono::high_resolution_clock::now();
    auto duration_rk = chrono::duration_cast<chrono::microseconds>(end_rk - start_rk);
    
    cout << "  Found " << rk_result.count << " occurrence(s) at positions: ";
    for (size_t i = 0; i < rk_result.positions.size(); ++i) {
        cout << rk_result.positions[i];
        if (i < rk_result.positions.size() - 1) cout << ", ";
    }
    cout << endl;
    cout << "  Time: " << duration_rk.count() << " microseconds" << endl << endl;
    
    // KMP
    cout << "Knuth-Morris-Pratt (KMP) Algorithm:" << endl;
    auto start_kmp = chrono::high_resolution_clock::now();
    SearchResult kmp_result = ClassicStringMatching::KMP::search(sequence, pattern);
    auto end_kmp = chrono::high_resolution_clock::now();
    auto duration_kmp = chrono::duration_cast<chrono::microseconds>(end_kmp - start_kmp);
    
    cout << "  Found " << kmp_result.count << " occurrence(s) at positions: ";
    for (size_t i = 0; i < kmp_result.positions.size(); ++i) {
        cout << kmp_result.positions[i];
        if (i < kmp_result.positions.size() - 1) cout << ", ";
    }
    cout << endl;
    cout << "  Time: " << duration_kmp.count() << " microseconds" << endl << endl;
    
    // Boyer-Moore
    cout << "Boyer-Moore Algorithm:" << endl;
    auto start_bm = chrono::high_resolution_clock::now();
    SearchResult bm_result = ClassicStringMatching::BoyerMoore::search(sequence, pattern);
    auto end_bm = chrono::high_resolution_clock::now();
    auto duration_bm = chrono::duration_cast<chrono::microseconds>(end_bm - start_bm);
    
    cout << "  Found " << bm_result.count << " occurrence(s) at positions: ";
    for (size_t i = 0; i < bm_result.positions.size(); ++i) {
        cout << bm_result.positions[i];
        if (i < bm_result.positions.size() - 1) cout << ", ";
    }
    cout << endl;
    cout << "  Time: " << duration_bm.count() << " microseconds" << endl << endl;
    
    // Comparison
    cout << "Algorithm Comparison:" << endl;
    cout << "  Rabin-Karp: " << duration_rk.count() << " μs" << endl;
    cout << "  KMP:        " << duration_kmp.count() << " μs" << endl;
    cout << "  Boyer-Moore: " << duration_bm.count() << " μs" << endl;
    
    // Verify all algorithms found same results
    bool results_match = (rk_result.count == kmp_result.count) && 
                        (kmp_result.count == bm_result.count);
    cout << "  Results consistent: " << (results_match ? "Yes" : "No") << endl << endl;
}

void demonstrateAdvancedAlgorithms() {
    printSeparator();
    cout << "ADVANCED ALGORITHMS DEMONSTRATION" << endl;
    printSeparator();
    
    string seq1 = "ATCGATCG";
    string seq2 = "ATCGATCC";
    
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl << endl;
    
    // Edit Distance Algorithms
    cout << "Edit Distance Algorithms:" << endl;
    cout << "  Levenshtein: " << EditDistanceAlgorithms::levenshteinDistance(seq1, seq2) << endl;
    cout << "  Damerau-Levenshtein: " << EditDistanceAlgorithms::damerauLevenshteinDistance(seq1, seq2) << endl;
    cout << "  DNA-specific: " << EditDistanceAlgorithms::dnaEditDistance(seq1, seq2) << endl;
    cout << "  Jaro-Winkler: " << fixed << setprecision(3) 
         << EditDistanceAlgorithms::jaroWinklerDistance(seq1, seq2) << endl;
    cout << "  LCS length: " << EditDistanceAlgorithms::longestCommonSubsequence(seq1, seq2) << endl << endl;
    
    // Embedding Search
    cout << "Embedding-Based Search:" << endl;
    EmbeddingSearch emb_search(64);
    emb_search.addSequence(seq1, 0);
    emb_search.addSequence(seq2, 1);
    emb_search.buildIndex();
    
    vector<pair<size_t, double>> emb_results = emb_search.search(seq1, 2);
    cout << "  Found " << emb_results.size() << " similar sequence(s)" << endl;
    for (const auto& result : emb_results) {
        cout << "    ID: " << result.first << ", Similarity: " 
             << fixed << setprecision(3) << result.second << endl;
    }
    cout << endl;
    
    // CNN Sequence Model
    cout << "CNN-Based Pattern Search:" << endl;
    CNNSequenceModel cnn_model(100, 32, 3);
    vector<size_t> cnn_positions = cnn_model.searchPattern(seq1 + seq1, "ATCG", 0.3);
    cout << "  Found pattern at " << cnn_positions.size() << " position(s)" << endl;
    if (!cnn_positions.empty()) {
        cout << "  Positions: ";
        for (size_t i = 0; i < cnn_positions.size() && i < 5; ++i) {
            cout << cnn_positions[i];
            if (i < cnn_positions.size() - 1 && i < 4) cout << ", ";
        }
        cout << endl;
    }
    cout << endl;
    
    // Lightweight LLM
    cout << "Lightweight LLM Search:" << endl;
    LightweightLLM llm(4, 64, 4, 2);
    vector<pair<size_t, double>> llm_results = llm.searchPattern(seq1 + seq1, "ATCG", 0.5);
    cout << "  Found " << llm_results.size() << " match(es)" << endl;
    for (size_t i = 0; i < llm_results.size() && i < 3; ++i) {
        cout << "    Position: " << llm_results[i].first 
             << ", Similarity: " << fixed << setprecision(3) << llm_results[i].second << endl;
    }
    cout << endl;
    
    // PIM Search
    cout << "PIM-Optimized Search:" << endl;
    PIMSearch pim_search;
    vector<string> sequences = {seq1 + seq1};
    pim_search.buildIndex(sequences);
    
    vector<size_t> pim_positions = pim_search.vectorizedSearch(seq1 + seq1, "ATCG");
    cout << "  Vectorized search found " << pim_positions.size() << " position(s)" << endl;
    
    vector<size_t> cache_positions = pim_search.cacheOptimizedSearch("ATCG");
    cout << "  Cache-optimized search found " << cache_positions.size() << " position(s)" << endl;
    cout << endl;
}

void demonstrateMCMCAndParallel() {
    printSeparator();
    cout << "MCMC AND PARALLEL SEARCH DEMONSTRATION" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    string initial_pattern = "ATCC";  // Slightly different from "ATCG"
    
    cout << "Sequence: " << sequence << endl;
    cout << "Initial Pattern: " << initial_pattern << endl << endl;
    
    // MCMC Search
    cout << "MCMC Pattern Evolution:" << endl;
    MCMCSearch mcmc(1000, 0.1, 1.0);
    MCMCSearch::MCMCResult mcmc_result = mcmc.search(sequence, initial_pattern);
    
    cout << "  Evolved Pattern: " << mcmc_result.evolved_pattern << endl;
    cout << "  Iterations: " << mcmc_result.iterations << endl;
    cout << "  Final Fitness: " << mcmc_result.final_fitness << endl;
    cout << "  Matches Found: " << mcmc_result.count << endl;
    if (mcmc_result.count > 0) {
        cout << "  Positions: ";
        for (size_t i = 0; i < mcmc_result.positions.size() && i < 5; ++i) {
            cout << mcmc_result.positions[i];
            if (i < mcmc_result.positions.size() - 1 && i < 4) cout << ", ";
        }
        cout << endl;
    }
    cout << endl;
    
    // Simulated Annealing
    cout << "Simulated Annealing Search:" << endl;
    MCMCSearch::MCMCResult sa_result = mcmc.simulatedAnnealingSearch(sequence, initial_pattern);
    cout << "  Evolved Pattern: " << sa_result.evolved_pattern << endl;
    cout << "  Matches Found: " << sa_result.count << endl << endl;
    
    // DDMCMC Search
    cout << "DDMCMC (Data-Driven MCMC) Search:" << endl;
    DDMCMCSearch ddmcmc(64, 500, 0.1, 1.0);
    vector<string> data_sequences = {sequence, "ATCGATCG", "GCTAGCTA", "TAGCTAGC"};
    ddmcmc.buildProposalDistribution(data_sequences);
    
    DDMCMCSearch::DDMCMCResult ddmcmc_result = ddmcmc.search("ATCG", 3);
    cout << "  Acceptance Rate: " << fixed << setprecision(3) 
         << ddmcmc_result.acceptance_rate << endl;
    cout << "  Matched Sequences: " << ddmcmc_result.matched_ids.size() << endl;
    for (size_t i = 0; i < ddmcmc_result.matched_ids.size(); ++i) {
        cout << "    ID: " << ddmcmc_result.matched_ids[i] 
             << ", Similarity: " << fixed << setprecision(3) 
             << ddmcmc_result.similarities[i] << endl;
    }
    cout << endl;
    
    // Parallel Search
    cout << "Parallel Search Methods:" << endl;
    string long_sequence = sequence + sequence + sequence + sequence;
    string pattern = "ATCG";
    
    auto start_par = chrono::high_resolution_clock::now();
    SearchResult par_result = ParallelSearch::parallelSearch(long_sequence, pattern, 4);
    auto end_par = chrono::high_resolution_clock::now();
    auto duration_par = chrono::duration_cast<chrono::microseconds>(end_par - start_par);
    
    cout << "  Parallel Search: " << par_result.count << " matches in " 
         << duration_par.count() << " μs" << endl;
    
    // Distributed Search
    auto start_dist = chrono::high_resolution_clock::now();
    SearchResult dist_result = ParallelSearch::distributedSearch(long_sequence, pattern, 500);
    auto end_dist = chrono::high_resolution_clock::now();
    auto duration_dist = chrono::duration_cast<chrono::microseconds>(end_dist - start_dist);
    
    cout << "  Distributed Search: " << dist_result.count << " matches in " 
         << duration_dist.count() << " μs" << endl;
    
    // Map-Reduce Search
    auto start_mr = chrono::high_resolution_clock::now();
    SearchResult mr_result = ParallelSearch::mapReduceSearch(long_sequence, pattern, 4);
    auto end_mr = chrono::high_resolution_clock::now();
    auto duration_mr = chrono::duration_cast<chrono::microseconds>(end_mr - start_mr);
    
    cout << "  Map-Reduce Search: " << mr_result.count << " matches in " 
         << duration_mr.count() << " μs" << endl;
    
    // Work-Stealing Search
    auto start_ws = chrono::high_resolution_clock::now();
    SearchResult ws_result = ParallelSearch::workStealingSearch(long_sequence, pattern, 4);
    auto end_ws = chrono::high_resolution_clock::now();
    auto duration_ws = chrono::duration_cast<chrono::microseconds>(end_ws - start_ws);
    
    cout << "  Work-Stealing Search: " << ws_result.count << " matches in " 
         << duration_ws.count() << " μs" << endl;
    cout << endl;
}

void demonstrateWARPCTC() {
    printSeparator();
    cout << "WARP-CTC PATTERN MATCHING DEMONSTRATION" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCG";
    string pattern = "ATCG";
    
    cout << "Sequence: " << sequence << endl;
    cout << "Pattern:  " << pattern << endl << endl;
    
    WarpCTC ctc(0.1);  // 10% blank probability
    
    // CTC Search
    cout << "WARP-CTC Search:" << endl;
    auto start_ctc = chrono::high_resolution_clock::now();
    WarpCTC::CTCResult ctc_result = ctc.search(sequence, pattern, 2);
    auto end_ctc = chrono::high_resolution_clock::now();
    auto duration_ctc = chrono::duration_cast<chrono::microseconds>(end_ctc - start_ctc);
    
    cout << "  Found " << ctc_result.num_matches << " alignment(s)" << endl;
    cout << "  Total Score: " << fixed << setprecision(4) << ctc_result.total_score << endl;
    if (!ctc_result.aligned_pattern.empty()) {
        cout << "  Best Alignment:" << endl;
        cout << "    Pattern: " << ctc_result.aligned_pattern << endl;
        cout << "    Sequence: " << ctc_result.aligned_sequence << endl;
    }
    cout << "  Time: " << duration_ctc.count() << " microseconds" << endl << endl;
    
    // Viterbi Decoding
    cout << "Viterbi Decoding (Best Path):" << endl;
    string viterbi_result = ctc.viterbiDecode(sequence, pattern);
    cout << "  Decoded Pattern: " << viterbi_result << endl << endl;
    
    // Beam Search
    cout << "Beam Search (Top Paths):" << endl;
    vector<pair<string, double>> beam_results = ctc.beamSearch(sequence, pattern, 3);
    cout << "  Found " << beam_results.size() << " path(s)" << endl;
    for (size_t i = 0; i < beam_results.size() && i < 3; ++i) {
        cout << "    Path " << (i+1) << ": " << beam_results[i].first 
             << " (score: " << fixed << setprecision(4) << beam_results[i].second << ")" << endl;
    }
    cout << endl;
    
    // Forward-Backward
    cout << "Forward-Backward Algorithm:" << endl;
    double forward_prob = ctc.forwardAlgorithm(sequence, pattern);
    double backward_prob = ctc.backwardAlgorithm(sequence, pattern);
    cout << "  Forward Probability: " << fixed << setprecision(6) << forward_prob << endl;
    cout << "  Backward Probability: " << fixed << setprecision(6) << backward_prob << endl;
    cout << "  Consistency Check: " << (std::abs(forward_prob - backward_prob) < 0.01 ? "Pass" : "Fail") << endl;
    cout << endl;
    
    // Alignment
    cout << "CTC Alignment:" << endl;
    auto aligned = ctc.align(sequence, pattern);
    cout << "  Aligned Pattern: " << aligned.first << endl;
    cout << "  Aligned Sequence: " << aligned.second << endl;
    cout << endl;
}

void demonstrateAdvancedAlgorithms2() {
    printSeparator();
    cout << "ADVANCED ALGORITHMS PART 2" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCG";
    vector<string> patterns = {"ATCG", "GCTA", "CGAT"};
    
    cout << "Sequence: " << sequence << endl;
    cout << "Patterns: ";
    for (const auto& p : patterns) cout << p << " ";
    cout << endl << endl;
    
    // Aho-Corasick
    cout << "Aho-Corasick Multi-Pattern Search:" << endl;
    AhoCorasick ac;
    ac.buildAutomaton(patterns);
    auto start_ac = chrono::high_resolution_clock::now();
    AhoCorasick::MultiPatternResult ac_result = ac.search(sequence);
    auto end_ac = chrono::high_resolution_clock::now();
    auto duration_ac = chrono::duration_cast<chrono::microseconds>(end_ac - start_ac);
    
    cout << "  Total Matches: " << ac_result.total_matches << endl;
    for (const auto& match_pair : ac_result.matches) {
        cout << "    " << match_pair.first << ": " << match_pair.second.size() << " occurrence(s)" << endl;
    }
    cout << "  Time: " << duration_ac.count() << " microseconds" << endl << endl;
    
    // Suffix Tree
    cout << "Suffix Tree Search:" << endl;
    SuffixTree st;
    st.build(sequence);
    auto start_st = chrono::high_resolution_clock::now();
    SearchResult st_result = st.search("ATCG");
    auto end_st = chrono::high_resolution_clock::now();
    auto duration_st = chrono::duration_cast<chrono::microseconds>(end_st - start_st);
    
    cout << "  Found " << st_result.count << " occurrence(s)" << endl;
    cout << "  Time: " << duration_st.count() << " microseconds" << endl << endl;
    
    // Suffix Array
    cout << "Suffix Array Search:" << endl;
    SuffixArray sa;
    sa.build(sequence);
    auto start_sa = chrono::high_resolution_clock::now();
    SearchResult sa_result = sa.search("ATCG");
    auto end_sa = chrono::high_resolution_clock::now();
    auto duration_sa = chrono::duration_cast<chrono::microseconds>(end_sa - start_sa);
    
    cout << "  Found " << sa_result.count << " occurrence(s)" << endl;
    string lrs = sa.longestRepeatedSubstring();
    cout << "  Longest Repeated Substring: " << lrs << endl;
    cout << "  Time: " << duration_sa.count() << " microseconds" << endl << endl;
    
    // Advanced DP
    cout << "Advanced Dynamic Programming:" << endl;
    string seq1 = "ACGTACGT";
    string seq2 = "TACGTACG";
    
    auto start_dp = chrono::high_resolution_clock::now();
    AlignmentResult dp_result = AdvancedDP::affineGapAlignment(seq1, seq2);
    auto end_dp = chrono::high_resolution_clock::now();
    auto duration_dp = chrono::duration_cast<chrono::microseconds>(end_dp - start_dp);
    
    cout << "  Affine Gap Alignment Score: " << dp_result.score << endl;
    cout << "  Time: " << duration_dp.count() << " microseconds" << endl << endl;
    
    // Advanced Deep Learning
    cout << "Advanced Deep Learning Approaches:" << endl;
    
    // LSTM
    DeepLearningAdvanced::LSTM lstm(4, 32, 1);
    vector<pair<size_t, double>> lstm_results = lstm.searchPattern(sequence, "ATCG", 0.3);
    cout << "  LSTM: Found " << lstm_results.size() << " match(es)" << endl;
    
    // GRU
    DeepLearningAdvanced::GRU gru(4, 32);
    vector<pair<size_t, double>> gru_results = gru.searchPattern(sequence, "ATCG", 0.3);
    cout << "  GRU: Found " << gru_results.size() << " match(es)" << endl;
    
    // Attention
    DeepLearningAdvanced::AttentionModel attention(64, 4);
    double attention_score = attention.alignSequences(seq1, seq2);
    cout << "  Attention Alignment Score: " << fixed << setprecision(3) << attention_score << endl;
    
    // Siamese Network
    DeepLearningAdvanced::SiameseNetwork siamese(64);
    vector<string> db_sequences = {sequence, "ATCGATCG", "GCTAGCTA"};
    vector<pair<size_t, double>> siamese_results = siamese.findSimilar("ATCG", db_sequences, 2);
    cout << "  Siamese Network: Found " << siamese_results.size() << " similar sequence(s)" << endl;
    for (const auto& result_pair : siamese_results) {
        cout << "    ID: " << result_pair.first << ", Similarity: " << fixed << setprecision(3) << result_pair.second << endl;
    }
    cout << endl;
}

void demonstrateHeuristicAndGraph() {
    printSeparator();
    cout << "HEURISTIC SEARCH AND GRAPH-BASED ALGORITHMS" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCG";
    vector<string> patterns = {"ATCG", "GCTA"};
    
    cout << "Sequence: " << sequence << endl;
    cout << "Patterns: ";
    for (const auto& p : patterns) cout << p << " ";
    cout << endl << endl;
    
    // Wu-Manber
    cout << "Wu-Manber Multi-Pattern Search:" << endl;
    WuManber wm(2);
    wm.buildTables(patterns);
    auto start_wm = chrono::high_resolution_clock::now();
    WuManber::MultiPatternResult wm_result = wm.search(sequence);
    auto end_wm = chrono::high_resolution_clock::now();
    auto duration_wm = chrono::duration_cast<chrono::microseconds>(end_wm - start_wm);
    
    cout << "  Total Matches: " << wm_result.total_matches << endl;
    for (const auto& match_pair : wm_result.matches) {
        cout << "    " << match_pair.first << ": " << match_pair.second.size() << " occurrence(s)" << endl;
    }
    cout << "  Time: " << duration_wm.count() << " microseconds" << endl << endl;
    
    // Bitap
    cout << "Bitap (Shift-Or) Algorithm:" << endl;
    Bitap bitap(1);  // Max 1 error
    auto start_bitap = chrono::high_resolution_clock::now();
    SearchResult bitap_result = bitap.search(sequence, "ATCC", 1);  // "ATCC" with 1 error matches "ATCG"
    auto end_bitap = chrono::high_resolution_clock::now();
    auto duration_bitap = chrono::duration_cast<chrono::microseconds>(end_bitap - start_bitap);
    
    cout << "  Found " << bitap_result.count << " match(es) with errors" << endl;
    cout << "  Time: " << duration_bitap.count() << " microseconds" << endl << endl;
    
    // Graph-based approaches
    cout << "Graph-Based Approaches:" << endl;
    
    // De Bruijn Graph
    GraphBased::DeBruijnGraph dbg(3);
    vector<string> sequences = {sequence, "ATCGATCG"};
    dbg.build(sequences);
    bool found = dbg.search("ATCG");
    cout << "  De Bruijn Graph: Pattern found = " << (found ? "Yes" : "No") << endl;
    string eulerian = dbg.findEulerianPath();
    cout << "  Eulerian Path: " << (eulerian.empty() ? "N/A" : eulerian.substr(0, 20)) << "..." << endl;
    
    // Overlap Graph
    GraphBased::OverlapGraph og;
    og.build(sequences, 2);
    string superstring = og.findShortestSuperstring();
    cout << "  Overlap Graph Superstring: " << superstring << endl;
    cout << endl;
    
    // Heuristic search
    cout << "Heuristic Search Algorithms:" << endl;
    
    string seq1 = "ACGTACGT";
    string seq2 = "TACGTACG";
    
    // A* Search
    auto start_astar = chrono::high_resolution_clock::now();
    int astar_score = HeuristicSearch::AStar::align(seq1, seq2);
    auto end_astar = chrono::high_resolution_clock::now();
    auto duration_astar = chrono::duration_cast<chrono::microseconds>(end_astar - start_astar);
    cout << "  A* Alignment Score: " << astar_score << endl;
    cout << "  Time: " << duration_astar.count() << " microseconds" << endl;
    
    // Beam Search
    auto start_beam = chrono::high_resolution_clock::now();
    int beam_score = HeuristicSearch::BeamSearch::align(seq1, seq2, 10);
    auto end_beam = chrono::high_resolution_clock::now();
    auto duration_beam = chrono::duration_cast<chrono::microseconds>(end_beam - start_beam);
    cout << "  Beam Search Score: " << beam_score << endl;
    cout << "  Time: " << duration_beam.count() << " microseconds" << endl;
    
    // Genetic Algorithm
    cout << "  Genetic Algorithm Pattern Evolution:" << endl;
    string evolved = HeuristicSearch::GeneticAlgorithm::evolvePattern(sequence, "ATCC", 30, 50, 0.1);
    cout << "    Evolved from 'ATCC' to: " << evolved << endl;
    
    // Simulated Annealing
    auto start_sa = chrono::high_resolution_clock::now();
    int sa_score = HeuristicSearch::SimulatedAnnealing::align(seq1, seq2, 100.0, 0.95, 500);
    auto end_sa = chrono::high_resolution_clock::now();
    auto duration_sa = chrono::duration_cast<chrono::microseconds>(end_sa - start_sa);
    cout << "  Simulated Annealing Score: " << sa_score << endl;
    cout << "  Time: " << duration_sa.count() << " microseconds" << endl;
    cout << endl;
}

void demonstrateVectorSimilarity() {
    printSeparator();
    cout << "VECTOR SIMILARITY ALGORITHMS" << endl;
    printSeparator();
    
    vector<string> sequences = {"ATCGATCG", "ATCGATCC", "GCTAGCTA", "TAGCTAGC", "ATCGATCGATCG"};
    string query = "ATCGATCG";
    
    cout << "Query Sequence: " << query << endl;
    cout << "Database Sequences: " << sequences.size() << endl << endl;
    
    // kNN
    cout << "k-Nearest Neighbors (kNN):" << endl;
    VectorSimilarity::KNN knn(3);
    knn.buildIndex(sequences);
    auto start_knn = chrono::high_resolution_clock::now();
    vector<VectorSimilarity::KNN::Neighbor> knn_results = knn.findNeighbors(query, 3);
    auto end_knn = chrono::high_resolution_clock::now();
    auto duration_knn = chrono::duration_cast<chrono::microseconds>(end_knn - start_knn);
    
    cout << "  Found " << knn_results.size() << " nearest neighbor(s):" << endl;
    for (const auto& neighbor : knn_results) {
        cout << "    Index: " << neighbor.index << ", Distance: " 
             << fixed << setprecision(4) << neighbor.distance 
             << ", Sequence: " << neighbor.sequence << endl;
    }
    cout << "  Time: " << duration_knn.count() << " microseconds" << endl << endl;
    
    // Approximate Nearest Neighbors
    cout << "Approximate Nearest Neighbors (ANN):" << endl;
    VectorSimilarity::ApproximateNearestNeighbors ann(10, 8);
    ann.buildIndex(sequences);
    auto start_ann = chrono::high_resolution_clock::now();
    vector<pair<size_t, double>> ann_results = ann.findNeighbors(query, 3);
    auto end_ann = chrono::high_resolution_clock::now();
    auto duration_ann = chrono::duration_cast<chrono::microseconds>(end_ann - start_ann);
    
    cout << "  Found " << ann_results.size() << " approximate neighbor(s):" << endl;
    for (const auto& result_pair : ann_results) {
        cout << "    Index: " << result_pair.first << ", Distance: " 
             << fixed << setprecision(4) << result_pair.second << endl;
    }
    cout << "  Time: " << duration_ann.count() << " microseconds" << endl << endl;
    
    // Dynamic Time Warping
    cout << "Dynamic Time Warping (DTW):" << endl;
    string seq1 = "ATCGATCG";
    string seq2 = "ATCGATCC";
    
    auto start_dtw = chrono::high_resolution_clock::now();
    double dtw_distance = VectorSimilarity::DynamicTimeWarping::distance(seq1, seq2);
    auto end_dtw = chrono::high_resolution_clock::now();
    auto duration_dtw = chrono::duration_cast<chrono::microseconds>(end_dtw - start_dtw);
    
    cout << "  DTW Distance: " << fixed << setprecision(4) << dtw_distance << endl;
    vector<pair<int, int>> dtw_path = VectorSimilarity::DynamicTimeWarping::findWarpingPath(seq1, seq2);
    cout << "  Warping Path Length: " << dtw_path.size() << " points" << endl;
    cout << "  Time: " << duration_dtw.count() << " microseconds" << endl << endl;
    
    // Vector Edit Distance
    cout << "Edit Distance in Vector Space:" << endl;
    auto start_ved = chrono::high_resolution_clock::now();
    double ved_distance = VectorSimilarity::VectorEditDistance::distance(seq1, seq2);
    auto end_ved = chrono::high_resolution_clock::now();
    auto duration_ved = chrono::duration_cast<chrono::microseconds>(end_ved - start_ved);
    
    cout << "  Vector Edit Distance: " << fixed << setprecision(4) << ved_distance << endl;
    cout << "  Time: " << duration_ved.count() << " microseconds" << endl << endl;
    
    // Vector Levenshtein Distance
    cout << "Levenshtein Distance in Vector Space:" << endl;
    auto start_vl = chrono::high_resolution_clock::now();
    double vl_distance = VectorSimilarity::VectorLevenshtein::distance(seq1, seq2);
    auto end_vl = chrono::high_resolution_clock::now();
    auto duration_vl = chrono::duration_cast<chrono::microseconds>(end_vl - start_vl);
    
    cout << "  Vector Levenshtein Distance: " << fixed << setprecision(4) << vl_distance << endl;
    cout << "  Time: " << duration_vl.count() << " microseconds" << endl << endl;
}

void demonstrateVariantCalling() {
    printSeparator();
    cout << "VARIANT CALLING" << endl;
    printSeparator();
    
    string reference = "ATCGATCGATCGATCG";
    string sequence = "ATCGATCCATCGATCG";  // SNP at position 7
    
    cout << "Reference: " << reference << endl;
    cout << "Sequence:  " << sequence << endl << endl;
    
    VariantCalling vc;
    
    // Single sequence variant calling
    cout << "Variant Calling (Single Sequence):" << endl;
    auto start_vc = chrono::high_resolution_clock::now();
    VariantCalling::VariantResult result = vc.callVariants(reference, sequence, 20.0);
    auto end_vc = chrono::high_resolution_clock::now();
    auto duration_vc = chrono::duration_cast<chrono::microseconds>(end_vc - start_vc);
    
    cout << "  Total Variants: " << result.total_variants << endl;
    cout << "  SNPs: " << result.snp_count << endl;
    cout << "  Indels: " << result.indel_count << endl;
    cout << "  Average Quality: " << fixed << setprecision(2) << result.average_quality << endl;
    
    for (const auto& variant : result.variants) {
        string type_str;
        switch (variant.type) {
            case VariantCalling::SNP: type_str = "SNP"; break;
            case VariantCalling::INSERTION: type_str = "INSERTION"; break;
            case VariantCalling::DELETION: type_str = "DELETION"; break;
            case VariantCalling::SUBSTITUTION: type_str = "SUBSTITUTION"; break;
            default: type_str = "COMPLEX"; break;
        }
        
        cout << "    " << type_str << " at position " << variant.position 
             << ": " << variant.ref_allele << " -> " << variant.alt_allele
             << " (Quality: " << fixed << setprecision(1) << variant.quality_score 
             << ", AF: " << fixed << setprecision(3) << variant.allele_frequency << ")" << endl;
        
        string annotation = VariantCalling::annotateVariant(variant, reference);
        cout << "      Annotation: " << annotation << endl;
    }
    cout << "  Time: " << duration_vc.count() << " microseconds" << endl << endl;
    
    // Pileup variant calling
    cout << "Variant Calling (Pileup):" << endl;
    vector<string> reads = {
        "ATCGATCCATCGATCG",
        "ATCGATCGATCGATCG",
        "ATCGATCCATCGATCG",
        "ATCGATCGATCGATCG",
        "ATCGATCCATCGATCG"
    };
    
    auto start_pileup = chrono::high_resolution_clock::now();
    VariantCalling::VariantResult pileup_result = vc.callVariantsPileup(reference, reads, 20.0, 3);
    auto end_pileup = chrono::high_resolution_clock::now();
    auto duration_pileup = chrono::duration_cast<chrono::microseconds>(end_pileup - start_pileup);
    
    cout << "  Reads: " << reads.size() << endl;
    cout << "  Variants Found: " << pileup_result.total_variants << endl;
    for (const auto& variant : pileup_result.variants) {
        cout << "    Position " << variant.position << ": " 
             << variant.ref_allele << " -> " << variant.alt_allele
             << " (Depth: " << variant.depth 
             << ", AF: " << fixed << setprecision(3) << variant.allele_frequency << ")" << endl;
    }
    cout << "  Time: " << duration_pileup.count() << " microseconds" << endl << endl;
}

void demonstrateMetagenomics() {
    printSeparator();
    cout << "METAGENOMICS ANALYSIS" << endl;
    printSeparator();
    
    // Simulate metagenomic reads from different organisms
    vector<string> reads = {
        "ATCGATCGATCG",  // Organism 1
        "ATCGATCGATCG",
        "GCTAGCTAGCTA",  // Organism 2
        "GCTAGCTAGCTA",
        "GCTAGCTAGCTA",
        "TAGCTAGCTAGC",  // Organism 3
        "TAGCTAGCTAGC"
    };
    
    // Reference database
    map<string, string> reference_db = {
        {"Bacillus_subtilis", "ATCGATCGATCG"},
        {"Escherichia_coli", "GCTAGCTAGCTA"},
        {"Staphylococcus_aureus", "TAGCTAGCTAGC"}
    };
    
    cout << "Metagenomic Sample:" << endl;
    cout << "  Total Reads: " << reads.size() << endl;
    cout << "  Reference Organisms: " << reference_db.size() << endl << endl;
    
    Metagenomics metagenomics;
    
    // Classify reads
    cout << "Read Classification:" << endl;
    auto start_classify = chrono::high_resolution_clock::now();
    Metagenomics::MetagenomicResult result = metagenomics.classifyReads(reads, reference_db, 0.7);
    auto end_classify = chrono::high_resolution_clock::now();
    auto duration_classify = chrono::duration_cast<chrono::microseconds>(end_classify - start_classify);
    
    cout << "  Classified Reads: " << result.classified_reads << " / " << result.total_reads << endl;
    cout << "  Organisms Found: " << result.organisms.size() << endl << endl;
    
    cout << "Organism Abundances:" << endl;
    for (const auto& org : result.organisms) {
        cout << "  " << org.organism_name << ":" << endl;
        cout << "    Abundance: " << fixed << setprecision(3) << org.abundance * 100 << "%" << endl;
        cout << "    Read Count: " << org.read_count << endl;
        cout << "    Average Similarity: " << fixed << setprecision(3) << org.similarity << endl;
    }
    cout << endl;
    
    cout << "Diversity Metrics:" << endl;
    cout << "  Shannon Diversity Index: " << fixed << setprecision(4) 
         << result.diversity_index << endl;
    
    map<string, double> abundances;
    for (const auto& org : result.organisms) {
        abundances[org.organism_name] = org.abundance;
    }
    double simpson = Metagenomics::calculateSimpsonDiversity(abundances);
    cout << "  Simpson Diversity Index: " << fixed << setprecision(4) << simpson << endl;
    cout << "  Time: " << duration_classify.count() << " microseconds" << endl << endl;
    
    // Marker gene analysis
    cout << "Marker Gene Analysis:" << endl;
    map<string, string> marker_db = {
        {"16S_rRNA_Bacillus", "ATCGATCG"},
        {"16S_rRNA_Ecoli", "GCTAGCTA"},
        {"16S_rRNA_Staph", "TAGCTAGC"}
    };
    
    auto start_marker = chrono::high_resolution_clock::now();
    Metagenomics::MetagenomicResult marker_result = metagenomics.assignTaxonomy(reads, marker_db);
    auto end_marker = chrono::high_resolution_clock::now();
    auto duration_marker = chrono::duration_cast<chrono::microseconds>(end_marker - start_marker);
    
    cout << "  Reads with Markers: " << marker_result.classified_reads << endl;
    cout << "  Taxonomic Assignments:" << endl;
    for (const auto& tax_pair : marker_result.taxonomy_counts) {
        cout << "    " << tax_pair.first << ": " << tax_pair.second << " reads" << endl;
    }
    cout << "  Time: " << duration_marker.count() << " microseconds" << endl << endl;
}

void demonstrateChordHashing() {
    printSeparator();
    cout << "CHORD DISTRIBUTED HASH TABLE (DHT)" << endl;
    printSeparator();
    
    // Create Chord ring
    ChordHashing chord(32);  // 32-bit identifier space
    
    cout << "Building Chord Ring:" << endl;
    
    // Add nodes
    vector<string> node_addresses = {"node1", "node2", "node3", "node4", "node5"};
    vector<size_t> node_ids;
    
    for (const string& addr : node_addresses) {
        size_t node_id = chord.addNode(addr);
        node_ids.push_back(node_id);
        cout << "  Added node: " << addr << " (ID: " << node_id << ")" << endl;
    }
    cout << endl;
    
    // Store sequences in DHT
    cout << "Storing DNA Sequences in Chord DHT:" << endl;
    vector<string> sequences = {
        "ATCGATCGATCG",
        "GCTAGCTAGCTA",
        "TAGCTAGCTAGC",
        "CGATCGATCGAT"
    };
    
    for (const string& seq : sequences) {
        auto start_store = chrono::high_resolution_clock::now();
        ChordHashing::LookupResult store_result = chord.store(seq, "metadata_" + seq);
        auto end_store = chrono::high_resolution_clock::now();
        auto duration_store = chrono::duration_cast<chrono::microseconds>(end_store - start_store);
        
        cout << "  Sequence: " << seq << endl;
        cout << "    Hash: " << store_result.key_hash << endl;
        cout << "    Stored at node: " << store_result.responsible_node_addr 
             << " (ID: " << store_result.responsible_node_id << ")" << endl;
        cout << "    Hops: " << store_result.hops << endl;
        cout << "    Time: " << duration_store.count() << " microseconds" << endl;
    }
    cout << endl;
    
    // Lookup sequences
    cout << "Looking up Sequences:" << endl;
    for (const string& seq : sequences) {
        auto start_lookup = chrono::high_resolution_clock::now();
        ChordHashing::LookupResult lookup_result = chord.lookup(seq);
        auto end_lookup = chrono::high_resolution_clock::now();
        auto duration_lookup = chrono::duration_cast<chrono::microseconds>(end_lookup - start_lookup);
        
        cout << "  Sequence: " << seq << endl;
        cout << "    Found at node: " << lookup_result.responsible_node_addr << endl;
        cout << "    Value: " << (lookup_result.value.empty() ? "Not found" : lookup_result.value) << endl;
        cout << "    Hops: " << lookup_result.hops << endl;
        cout << "    Time: " << duration_lookup.count() << " microseconds" << endl;
    }
    cout << endl;
    
    // Display node data distribution
    cout << "Data Distribution:" << endl;
    for (size_t node_id : node_ids) {
        auto node_data = chord.getNodeData(node_id);
        cout << "  Node " << node_id << ": " << node_data.size() << " key(s)" << endl;
        for (const auto& data_pair : node_data) {
            cout << "    Key hash: " << data_pair.first << " -> " << data_pair.second << endl;
        }
    }
    cout << endl;
    
    // Ring statistics
    cout << "Chord Ring Statistics:" << endl;
    ChordHashing::RingStatistics stats = chord.getStatistics();
    cout << "  Number of Nodes: " << stats.num_nodes << endl;
    cout << "  Total Keys: " << stats.total_keys << endl;
    cout << "  Average Keys per Node: " << fixed << setprecision(2) 
         << stats.avg_keys_per_node << endl;
    cout << "  Max Keys per Node: " << stats.max_keys_per_node << endl;
    cout << "  Min Keys per Node: " << stats.min_keys_per_node << endl;
    cout << "  Average Finger Table Size: " << fixed << setprecision(2) 
         << stats.avg_finger_table_size << endl;
    cout << endl;
    
    // Demonstrate node removal
    cout << "Node Removal (Redistributing Data):" << endl;
    if (!node_ids.empty()) {
        size_t removed_node = node_ids[0];
        cout << "  Removing node: " << removed_node << endl;
        chord.removeNode(removed_node);
        
        ChordHashing::RingStatistics new_stats = chord.getStatistics();
        cout << "  Nodes after removal: " << new_stats.num_nodes << endl;
        cout << "  Keys redistributed: " << new_stats.total_keys << endl;
    }
    cout << endl;
}

void demonstrateProbabilisticML() {
    printSeparator();
    cout << "PROBABILISTIC MACHINE LEARNING METHODS" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCG";
    vector<string> patterns = {"ATCG", "GCTA", "TAGC"};
    
    // Naive Bayes
    cout << "Naive Bayes Classification:" << endl;
    ProbabilisticML::NaiveBayes nb;
    vector<string> training_seqs = {sequence, "GCTAGCTAGCTA", "TAGCTAGCTAGC"};
    vector<string> training_labels = {"pattern1", "pattern2", "pattern3"};
    nb.train(training_seqs, training_labels);
    
    auto start_nb = chrono::high_resolution_clock::now();
    ProbabilisticML::NaiveBayes::ClassificationResult nb_result = nb.classify("ATCGATCG");
    auto end_nb = chrono::high_resolution_clock::now();
    auto duration_nb = chrono::duration_cast<chrono::microseconds>(end_nb - start_nb);
    
    cout << "  Query: ATCGATCG" << endl;
    cout << "  Predicted Class: " << nb_result.predicted_class << endl;
    cout << "  Confidence: " << fixed << setprecision(3) << nb_result.confidence << endl;
    cout << "  Time: " << duration_nb.count() << " microseconds" << endl << endl;
    
    // Restricted Boltzmann Machine
    cout << "Restricted Boltzmann Machine (RBM):" << endl;
    ProbabilisticML::RestrictedBoltzmannMachine rbm(64, 32);
    
    auto start_rbm = chrono::high_resolution_clock::now();
    rbm.train(training_seqs);
    auto end_rbm = chrono::high_resolution_clock::now();
    auto duration_rbm = chrono::duration_cast<chrono::microseconds>(end_rbm - start_rbm);
    
    cout << "  Training Time: " << duration_rbm.count() << " microseconds" << endl;
    
    vector<double> rbm_features = rbm.encode(sequence);
    cout << "  Encoded Features: " << rbm_features.size() << " dimensions" << endl;
    cout << "  Sample Features: ";
    for (size_t i = 0; i < std::min(rbm_features.size(), size_t(5)); ++i) {
        cout << fixed << setprecision(3) << rbm_features[i] << " ";
    }
    cout << "..." << endl;
    
    auto start_sim = chrono::high_resolution_clock::now();
    vector<pair<size_t, double>> similar = rbm.findSimilar(sequence, training_seqs, 3);
    auto end_sim = chrono::high_resolution_clock::now();
    auto duration_sim = chrono::duration_cast<chrono::microseconds>(end_sim - start_sim);
    
    cout << "  Similar Sequences: " << similar.size() << endl;
    for (const auto& sim_pair : similar) {
        cout << "    Index " << sim_pair.first << ": Similarity " 
             << fixed << setprecision(3) << sim_pair.second << endl;
    }
    cout << "  Search Time: " << duration_sim.count() << " microseconds" << endl << endl;
    
    // Deep Belief Network
    cout << "Deep Belief Network (DBN):" << endl;
    vector<int> layer_sizes = {32, 16, 8};
    ProbabilisticML::DeepBeliefNetwork dbn(layer_sizes);
    
    auto start_dbn = chrono::high_resolution_clock::now();
    dbn.train(training_seqs);
    auto end_dbn = chrono::high_resolution_clock::now();
    auto duration_dbn = chrono::duration_cast<chrono::microseconds>(end_dbn - start_dbn);
    
    cout << "  Training Time: " << duration_dbn.count() << " microseconds" << endl;
    
    vector<double> dbn_features = dbn.encode(sequence);
    cout << "  Encoded Features: " << dbn_features.size() << " dimensions" << endl;
    
    auto start_match = chrono::high_resolution_clock::now();
    vector<pair<string, double>> dbn_matches = dbn.matchPatterns(sequence, patterns);
    auto end_match = chrono::high_resolution_clock::now();
    auto duration_match = chrono::duration_cast<chrono::microseconds>(end_match - start_match);
    
    cout << "  Pattern Matches:" << endl;
    for (const auto& match_pair : dbn_matches) {
        cout << "    " << match_pair.first << ": Similarity " 
             << fixed << setprecision(3) << match_pair.second << endl;
    }
    cout << "  Match Time: " << duration_match.count() << " microseconds" << endl << endl;
    
    // Markov Decision Process
    cout << "Markov Decision Process (MDP):" << endl;
    ProbabilisticML::MarkovDecisionProcess mdp(2.0, -1.0, -1.0);
    
    auto start_mdp = chrono::high_resolution_clock::now();
    ProbabilisticML::MarkovDecisionProcess::MDPResult mdp_result = mdp.findPattern(sequence, "ATCG", 1);
    auto end_mdp = chrono::high_resolution_clock::now();
    auto duration_mdp = chrono::duration_cast<chrono::microseconds>(end_mdp - start_mdp);
    
    cout << "  Pattern: ATCG" << endl;
    cout << "  Matched Positions: " << mdp_result.positions.size() << endl;
    for (int pos : mdp_result.positions) {
        cout << "    Position: " << pos << endl;
    }
    cout << "  Total Reward: " << fixed << setprecision(3) << mdp_result.total_reward << endl;
    cout << "  Actions Taken: " << mdp_result.action_sequence.size() << endl;
    cout << "  Time: " << duration_mdp.count() << " microseconds" << endl << endl;
    
    // Markov Random Field
    cout << "Markov Random Field (MRF):" << endl;
    ProbabilisticML::MarkovRandomField mrf(-2.0, 1.0, 0.5);
    
    auto start_mrf = chrono::high_resolution_clock::now();
    ProbabilisticML::MarkovRandomField::MRFResult mrf_result = mrf.findPattern(sequence, "ATCG");
    auto end_mrf = chrono::high_resolution_clock::now();
    auto duration_mrf = chrono::duration_cast<chrono::microseconds>(end_mrf - start_mrf);
    
    cout << "  Pattern: ATCG" << endl;
    cout << "  Matched Positions: " << mrf_result.matched_positions.size() << endl;
    for (int pos : mrf_result.matched_positions) {
        cout << "    Position: " << pos << endl;
    }
    cout << "  Energy: " << fixed << setprecision(3) << mrf_result.energy << endl;
    cout << "  Time: " << duration_mrf.count() << " microseconds" << endl;
    
    // Belief Propagation
    auto start_bp = chrono::high_resolution_clock::now();
    ProbabilisticML::MarkovRandomField::MRFResult bp_result = mrf.beliefPropagation(sequence, "ATCG");
    auto end_bp = chrono::high_resolution_clock::now();
    auto duration_bp = chrono::duration_cast<chrono::microseconds>(end_bp - start_bp);
    
    cout << "  Belief Propagation:" << endl;
    cout << "    Matched Positions: " << bp_result.matched_positions.size() << endl;
    cout << "    Energy: " << fixed << setprecision(3) << bp_result.energy << endl;
    cout << "    Time: " << duration_bp.count() << " microseconds" << endl << endl;
}

void demonstrateConcurrentMultiSearch() {
    printSeparator();
    cout << "CONCURRENT MULTI-TECHNIQUE SUBSTRING SEARCH" << endl;
    printSeparator();
    
    string sequence = "ATCGATCGATCGATCGATCGATCG";
    string pattern = "ATCG";
    
    // Default configuration (fast exact techniques)
    cout << "Default Configuration (Fast Exact Techniques):" << endl;
    ConcurrentMultiSearch::SearchConfig default_config = ConcurrentMultiSearch::getDefaultConfig();
    
    auto start_default = chrono::high_resolution_clock::now();
    ConcurrentMultiSearch searcher;
    ConcurrentMultiSearch::MultiTechniqueResult default_result = searcher.search(sequence, pattern, default_config);
    auto end_default = chrono::high_resolution_clock::now();
    auto duration_default = chrono::duration_cast<chrono::microseconds>(end_default - start_default);
    
    cout << "  Pattern: " << pattern << endl;
    cout << "  Techniques Used: " << default_result.num_techniques_used << endl;
    cout << "  Total Unique Matches: " << default_result.total_unique_matches << endl;
    cout << "  Consensus Positions: " << default_result.consensus_positions.size() << endl;
    cout << "  Total Time: " << duration_default.count() << " microseconds" << endl;
    
    cout << "  Results by Technique:" << endl;
    for (const auto& tech_pair : default_result.technique_counts) {
        cout << "    " << tech_pair.first << ": " << tech_pair.second 
             << " matches, " << fixed << setprecision(2) 
             << default_result.technique_times[tech_pair.first] << " μs" << endl;
    }
    cout << endl;
    
    // Comprehensive configuration (all techniques including fuzzy)
    cout << "Comprehensive Configuration (All Techniques):" << endl;
    ConcurrentMultiSearch::SearchConfig comp_config = ConcurrentMultiSearch::getComprehensiveConfig();
    
    auto start_comp = chrono::high_resolution_clock::now();
    ConcurrentMultiSearch::MultiTechniqueResult comp_result = searcher.search(sequence, pattern, comp_config);
    auto end_comp = chrono::high_resolution_clock::now();
    auto duration_comp = chrono::duration_cast<chrono::microseconds>(end_comp - start_comp);
    
    cout << "  Pattern: " << pattern << endl;
    cout << "  Techniques Used: " << comp_result.num_techniques_used << endl;
    cout << "  Total Unique Matches: " << comp_result.total_unique_matches << endl;
    cout << "  Consensus Positions (≥2 techniques): " << comp_result.consensus_positions.size() << endl;
    if (!comp_result.consensus_positions.empty()) {
        cout << "    Positions: ";
        for (size_t i = 0; i < std::min(comp_result.consensus_positions.size(), size_t(10)); ++i) {
            cout << comp_result.consensus_positions[i] << " ";
        }
        if (comp_result.consensus_positions.size() > 10) {
            cout << "...";
        }
        cout << endl;
    }
    cout << "  Unique Positions (single technique): " << comp_result.unique_positions.size() << endl;
    cout << "  Total Time: " << duration_comp.count() << " microseconds" << endl;
    
    cout << "  Results by Technique:" << endl;
    for (const auto& tech_pair : comp_result.technique_counts) {
        cout << "    " << tech_pair.first << ": " << tech_pair.second 
             << " matches, " << fixed << setprecision(2) 
             << comp_result.technique_times[tech_pair.first] << " μs" << endl;
    }
    cout << endl;
    
    // Custom configuration (specific techniques)
    cout << "Custom Configuration (KMP + Boyer-Moore):" << endl;
    ConcurrentMultiSearch::SearchConfig custom_config;
    custom_config.techniques = {ConcurrentMultiSearch::KMP, ConcurrentMultiSearch::BOYER_MOORE};
    custom_config.use_threads = true;
    custom_config.num_threads = 2;
    
    auto start_custom = chrono::high_resolution_clock::now();
    ConcurrentMultiSearch::MultiTechniqueResult custom_result = searcher.search(sequence, pattern, custom_config);
    auto end_custom = chrono::high_resolution_clock::now();
    auto duration_custom = chrono::duration_cast<chrono::microseconds>(end_custom - start_custom);
    
    cout << "  Pattern: " << pattern << endl;
    cout << "  Techniques Used: " << custom_result.num_techniques_used << endl;
    cout << "  Total Unique Matches: " << custom_result.total_unique_matches << endl;
    cout << "  Total Time: " << duration_custom.count() << " microseconds" << endl;
    
    cout << "  Results by Technique:" << endl;
    for (const auto& tech_pair : custom_result.technique_counts) {
        cout << "    " << tech_pair.first << ": " << tech_pair.second 
             << " matches, " << fixed << setprecision(2) 
             << custom_result.technique_times[tech_pair.first] << " μs" << endl;
    }
    cout << endl;
    
    // Performance comparison: sequential vs concurrent
    cout << "Performance Comparison (Sequential vs Concurrent):" << endl;
    
    // Sequential execution
    auto start_seq = chrono::high_resolution_clock::now();
    ExactMatch exact;
    NaiveSearch naive;
    ClassicStringMatching::KMP kmp;
    ClassicStringMatching::BoyerMoore bm;
    exact.search(sequence, pattern);
    naive.search(sequence, pattern);
    kmp.search(sequence, pattern);
    bm.search(sequence, pattern);
    auto end_seq = chrono::high_resolution_clock::now();
    auto duration_seq = chrono::duration_cast<chrono::microseconds>(end_seq - start_seq);
    
    cout << "  Sequential Time: " << duration_seq.count() << " microseconds" << endl;
    cout << "  Concurrent Time: " << duration_default.count() << " microseconds" << endl;
    cout << "  Speedup: " << fixed << setprecision(2) 
         << static_cast<double>(duration_seq.count()) / duration_default.count() << "x" << endl;
    cout << endl;
}

void demonstrateSkipGraph() {
    printSeparator();
    cout << "SKIP-GRAPH HIERARCHICAL INDEXING FOR DNA SEQUENCES" << endl;
    printSeparator();
    
    // Generate a longer sequence for demonstration
    string long_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    for (int i = 0; i < 5; ++i) {
        long_sequence += long_sequence;
    }
    
    cout << "Sequence Length: " << long_sequence.length() << " bases" << endl;
    cout << "Building Skip-Graph Index..." << endl;
    
    // Default configuration
    SkipGraph::GraphConfig config;
    config.max_levels = 16;
    config.subsequence_length = 4;
    config.use_hierarchical = true;
    config.hash_function = 1; // Rolling hash
    
    auto start_build = chrono::high_resolution_clock::now();
    SkipGraph skip_graph(config);
    skip_graph.buildIndex(long_sequence);
    auto end_build = chrono::high_resolution_clock::now();
    auto duration_build = chrono::duration_cast<chrono::microseconds>(end_build - start_build);
    
    cout << "Index Built in: " << duration_build.count() << " microseconds" << endl;
    
    // Get statistics
    SkipGraph::Statistics stats = skip_graph.getStatistics();
    cout << "\nSkip-Graph Statistics:" << endl;
    cout << "  Total Nodes: " << stats.total_nodes << endl;
    cout << "  Total Subsequences Indexed: " << stats.total_subsequences << endl;
    cout << "  Maximum Level: " << stats.max_level << endl;
    cout << "  Estimated Memory Usage: " << stats.memory_usage << " bytes" << endl;
    cout << "  Nodes per Level:" << endl;
    for (const auto& level_pair : stats.nodes_per_level) {
        cout << "    Level " << level_pair.first << ": " << level_pair.second << " nodes" << endl;
    }
    cout << endl;
    
    // Search for patterns
    vector<string> patterns = {"ATCG", "CGAT", "GATC", "TCGA"};
    
    cout << "Searching for Patterns:" << endl;
    for (const string& pattern : patterns) {
        SkipGraph::SkipGraphResult result = skip_graph.search(pattern);
        
        cout << "  Pattern: " << pattern << endl;
        cout << "    Found: " << (result.found ? "Yes" : "No") << endl;
        cout << "    Matches: " << result.num_matches << endl;
        cout << "    Levels Searched: " << result.levels_searched << endl;
        cout << "    Search Time: " << fixed << setprecision(2) << result.search_time << " μs" << endl;
        if (result.num_matches > 0 && result.num_matches <= 10) {
            cout << "    Positions: ";
            for (int pos : result.positions) {
                cout << pos << " ";
            }
            cout << endl;
        } else if (result.num_matches > 10) {
            cout << "    First 10 Positions: ";
            for (size_t i = 0; i < 10; ++i) {
                cout << result.positions[i] << " ";
            }
            cout << "..." << endl;
        }
        cout << endl;
    }
    
    // Hierarchical search demonstration
    cout << "Hierarchical Search (Level-by-Level):" << endl;
    string test_pattern = "ATCG";
    for (int level = 0; level <= stats.max_level && level < 5; ++level) {
        SkipGraph::SkipGraphResult result = skip_graph.searchAtLevel(test_pattern, level);
        cout << "  Level " << level << ": " 
             << (result.found ? "Found" : "Not Found") 
             << " (" << result.search_time << " μs)" << endl;
    }
    cout << endl;
    
    // Pre-caching demonstration
    cout << "Pre-caching Additional Patterns..." << endl;
    skip_graph.preCacheSubsequences("ATCGATCGATCG", 6);
    SkipGraph::Statistics new_stats = skip_graph.getStatistics();
    cout << "  New Total Subsequences: " << new_stats.total_subsequences << endl;
    cout << endl;
    
    // Performance comparison: indexed vs non-indexed
    cout << "Performance Comparison:" << endl;
    string search_pattern = "ATCG";
    
    // Indexed search
    auto start_indexed = chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; ++i) {
        skip_graph.search(search_pattern);
    }
    auto end_indexed = chrono::high_resolution_clock::now();
    auto duration_indexed = chrono::duration_cast<chrono::microseconds>(end_indexed - start_indexed);
    
    // Non-indexed search (naive)
    auto start_naive = chrono::high_resolution_clock::now();
    NaiveSearch naive;
    for (int i = 0; i < 100; ++i) {
        naive.search(long_sequence, search_pattern);
    }
    auto end_naive = chrono::high_resolution_clock::now();
    auto duration_naive = chrono::duration_cast<chrono::microseconds>(end_naive - start_naive);
    
    cout << "  100 Searches:" << endl;
    cout << "    Indexed (Skip-Graph): " << duration_indexed.count() << " μs" << endl;
    cout << "    Non-Indexed (Naive): " << duration_naive.count() << " μs" << endl;
    if (duration_naive.count() > 0) {
        cout << "    Speedup: " << fixed << setprecision(2) 
             << static_cast<double>(duration_naive.count()) / duration_indexed.count() << "x" << endl;
    }
    cout << endl;
}

void demonstrateDancingLinks() {
    printSeparator();
    cout << "DANCING LINKS (ALGORITHM X) FOR EXACT COVER ON SPARSE-ENTROPIC DNA" << endl;
    printSeparator();
    
    // Create a sparse-entropic (low entropy, highly repetitive) DNA sequence
    string sparse_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    for (int i = 0; i < 3; ++i) {
        sparse_sequence += sparse_sequence;
    }
    
    cout << "Sparse-Entropic Sequence Length: " << sparse_sequence.length() << " bases" << endl;
    cout << "Sequence: " << sparse_sequence.substr(0, 50) << "..." << endl;
    
    // Test 1: Exact cover with single pattern
    cout << "\nTest 1: Exact Cover with Single Pattern" << endl;
    string pattern = "ATCG";
    
    DancingLinks::DLXConfig config;
    config.find_all_solutions = false;
    config.max_solutions = 1;
    config.sparse_mode = true;
    config.low_entropy_mode = true;
    
    DancingLinks dlx(config);
    dlx.buildSparseEntropicProblem(sparse_sequence, pattern);
    
    DancingLinks::Statistics stats = dlx.getStatistics();
    cout << "  Problem Statistics:" << endl;
    cout << "    Rows (Pattern Occurrences): " << stats.num_rows << endl;
    cout << "    Columns (Sequence Positions): " << stats.num_columns << endl;
    cout << "    Nodes: " << stats.num_nodes << endl;
    cout << "    Density: " << fixed << setprecision(4) << stats.density << endl;
    cout << "    Entropy: " << stats.entropy << endl;
    
    auto start_solve = chrono::high_resolution_clock::now();
    DancingLinks::ExactCoverSolution solution = dlx.solve();
    auto end_solve = chrono::high_resolution_clock::now();
    auto duration_solve = chrono::duration_cast<chrono::microseconds>(end_solve - start_solve);
    
    cout << "  Solution:" << endl;
    cout << "    Found: " << (solution.found ? "Yes" : "No") << endl;
    cout << "    Selected Rows: " << solution.selected_rows.size() << endl;
    cout << "    Solve Time: " << solution.solve_time << " μs" << endl;
    if (solution.selected_rows.size() <= 10) {
        cout << "    Row Indices: ";
        for (int row : solution.selected_rows) {
            cout << row << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    // Test 2: Exact cover with multiple patterns
    cout << "Test 2: Exact Cover with Multiple Patterns" << endl;
    vector<string> patterns = {"ATCG", "CGAT", "GATC"};
    
    DancingLinks dlx2(config);
    dlx2.buildProblem(sparse_sequence, patterns);
    
    DancingLinks::Statistics stats2 = dlx2.getStatistics();
    cout << "  Problem Statistics:" << endl;
    cout << "    Rows: " << stats2.num_rows << endl;
    cout << "    Columns: " << stats2.num_columns << endl;
    cout << "    Nodes: " << stats2.num_nodes << endl;
    cout << "    Density: " << fixed << setprecision(4) << stats2.density << endl;
    cout << "    Entropy: " << stats2.entropy << endl;
    
    DancingLinks::ExactCoverSolution solution2 = dlx2.solve();
    cout << "  Solution:" << endl;
    cout << "    Found: " << (solution2.found ? "Yes" : "No") << endl;
    cout << "    Selected Rows: " << solution2.selected_rows.size() << endl;
    cout << "    Solve Time: " << solution2.solve_time << " μs" << endl;
    cout << endl;
    
    // Test 3: Find all solutions
    cout << "Test 3: Find All Solutions (Limited)" << endl;
    DancingLinks::DLXConfig config_all;
    config_all.find_all_solutions = true;
    config_all.max_solutions = 5;
    config_all.sparse_mode = true;
    config_all.low_entropy_mode = true;
    
    DancingLinks dlx3(config_all);
    string short_sequence = "ATCGATCGATCG";
    dlx3.buildSparseEntropicProblem(short_sequence, "ATCG");
    
    vector<DancingLinks::ExactCoverSolution> all_solutions = dlx3.getAllSolutions();
    cout << "  Total Solutions Found: " << all_solutions.size() << endl;
    for (size_t i = 0; i < all_solutions.size() && i < 3; ++i) {
        cout << "    Solution " << (i + 1) << ": " << all_solutions[i].selected_rows.size() << " rows" << endl;
    }
    cout << endl;
    
    // Test 4: Performance comparison
    cout << "Test 4: Performance on Different Sequence Types" << endl;
    
    // Low entropy sequence
    string low_entropy = "AAAAAAAATTTTTTTTCCCCCCCCGGGGGGGG";
    for (int i = 0; i < 2; ++i) {
        low_entropy += low_entropy;
    }
    
    // High entropy sequence
    string high_entropy = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    for (int i = 0; i < 2; ++i) {
        high_entropy += high_entropy;
    }
    
    DancingLinks dlx_low(config);
    dlx_low.buildSparseEntropicProblem(low_entropy, "AAAA");
    auto start_low = chrono::high_resolution_clock::now();
    DancingLinks::ExactCoverSolution sol_low = dlx_low.solve();
    auto end_low = chrono::high_resolution_clock::now();
    auto duration_low = chrono::duration_cast<chrono::microseconds>(end_low - start_low);
    
    DancingLinks dlx_high(config);
    dlx_high.buildSparseEntropicProblem(high_entropy, "ATCG");
    auto start_high = chrono::high_resolution_clock::now();
    DancingLinks::ExactCoverSolution sol_high = dlx_high.solve();
    auto end_high = chrono::high_resolution_clock::now();
    auto duration_high = chrono::duration_cast<chrono::microseconds>(end_high - start_high);
    
    cout << "  Low Entropy Sequence:" << endl;
    cout << "    Length: " << low_entropy.length() << endl;
    cout << "    Entropy: " << dlx_low.getStatistics().entropy << endl;
    cout << "    Solve Time: " << duration_low.count() << " μs" << endl;
    cout << "    Solution Found: " << (sol_low.found ? "Yes" : "No") << endl;
    
    cout << "  High Entropy Sequence:" << endl;
    cout << "    Length: " << high_entropy.length() << endl;
    cout << "    Entropy: " << dlx_high.getStatistics().entropy << endl;
    cout << "    Solve Time: " << duration_high.count() << " μs" << endl;
    cout << "    Solution Found: " << (sol_high.found ? "Yes" : "No") << endl;
    cout << endl;
}

void printMenu() {
    cout << "\nDNA Sequence Alignment Tool" << endl;
    cout << "============================" << endl;
    cout << "1. Exact Match Search" << endl;
    cout << "2. Naive Search" << endl;
    cout << "3. Smith-Waterman (Local Alignment)" << endl;
    cout << "4. Needleman-Wunsch (Global Alignment)" << endl;
    cout << "5. Fuzzy Search (Edit Distance)" << endl;
    cout << "6. Compression (Grammar & Lossy)" << endl;
    cout << "7. Classic String Matching (Rabin-Karp, KMP, Boyer-Moore)" << endl;
    cout << "8. Advanced Algorithms (Edit Dist, Embedding, CNN, LLM, PIM)" << endl;
    cout << "9. MCMC & Parallel Search (Pattern Evolution, DDMCMC, Distributed)" << endl;
    cout << "10. WARP-CTC Pattern Matching" << endl;
    cout << "11. Advanced Algorithms 2 (Aho-Corasick, Suffix Tree/Array, Advanced DP, DL)" << endl;
    cout << "12. Heuristic & Graph-Based (Wu-Manber, Bitap, Graph, A*, GA, SA)" << endl;
    cout << "13. Vector Similarity (kNN, ANN, DTW, Vector Edit/Levenshtein)" << endl;
    cout << "14. Variant Calling" << endl;
    cout << "15. Metagenomics Analysis" << endl;
    cout << "16. Chord Distributed Hash Table (DHT)" << endl;
    cout << "17. Probabilistic ML (Naive Bayes, DBN, RBM, MDP, MRF)" << endl;
    cout << "18. Concurrent Multi-Technique Search" << endl;
    cout << "19. Skip-Graph Hierarchical Indexing" << endl;
    cout << "20. Dancing Links (Algorithm X) for Exact Cover" << endl;
    cout << "21. Compare Local vs Global Alignment" << endl;
    cout << "22. Run All Demonstrations" << endl;
    cout << "0. Exit" << endl;
    cout << "\nSelect an option: ";
}

int main(int argc, char* argv[]) {
    // If command line arguments provided, use them
    if (argc > 1) {
        string mode = argv[1];
        
        if (mode == "exact") {
            demonstrateExactMatch();
        } else if (mode == "naive") {
            demonstrateNaiveSearch();
        } else if (mode == "smith") {
            demonstrateSmithWaterman();
        } else if (mode == "needleman") {
            demonstrateNeedlemanWunsch();
        } else if (mode == "fuzzy") {
            demonstrateFuzzySearch();
        } else if (mode == "compress") {
            demonstrateCompression();
        } else if (mode == "classic") {
            demonstrateClassicStringMatching();
        } else if (mode == "advanced") {
            demonstrateAdvancedAlgorithms();
        } else if (mode == "mcmc") {
            demonstrateMCMCAndParallel();
        } else if (mode == "ctc") {
            demonstrateWARPCTC();
        } else if (mode == "advanced2") {
            demonstrateAdvancedAlgorithms2();
        } else if (mode == "heuristic") {
            demonstrateHeuristicAndGraph();
        } else if (mode == "vector") {
            demonstrateVectorSimilarity();
        } else if (mode == "variant") {
            demonstrateVariantCalling();
        } else if (mode == "metagenomics") {
            demonstrateMetagenomics();
        } else if (mode == "chord") {
            demonstrateChordHashing();
        } else if (mode == "probabilistic" || mode == "ml") {
            demonstrateProbabilisticML();
        } else if (mode == "concurrent" || mode == "multi") {
            demonstrateConcurrentMultiSearch();
        } else if (mode == "skip" || mode == "skipgraph") {
            demonstrateSkipGraph();
        } else if (mode == "dancing" || mode == "dlx") {
            demonstrateDancingLinks();
        } else if (mode == "compare") {
            demonstrateComparison();
        } else if (mode == "all") {
            demonstrateExactMatch();
            demonstrateNaiveSearch();
            demonstrateSmithWaterman();
            demonstrateNeedlemanWunsch();
            demonstrateFuzzySearch();
            demonstrateCompression();
            demonstrateClassicStringMatching();
            demonstrateAdvancedAlgorithms();
            demonstrateMCMCAndParallel();
            demonstrateWARPCTC();
            demonstrateAdvancedAlgorithms2();
            demonstrateHeuristicAndGraph();
            demonstrateVectorSimilarity();
            demonstrateVariantCalling();
            demonstrateMetagenomics();
            demonstrateChordHashing();
            demonstrateProbabilisticML();
            demonstrateConcurrentMultiSearch();
            demonstrateSkipGraph();
            demonstrateDancingLinks();
            demonstrateComparison();
        } else {
            cout << "Usage: " << argv[0] << " [exact|naive|smith|needleman|fuzzy|compress|classic|advanced|mcmc|ctc|advanced2|heuristic|vector|variant|metagenomics|chord|probabilistic|ml|concurrent|multi|skip|skipgraph|dancing|dlx|compare|all]" << endl;
            cout << "Or run without arguments for interactive mode." << endl;
            return 1;
        }
        return 0;
    }
    
    // Interactive mode
    int choice;
    do {
        printMenu();
        cin >> choice;
        
        switch (choice) {
            case 1:
                demonstrateExactMatch();
                break;
            case 2:
                demonstrateNaiveSearch();
                break;
            case 3:
                demonstrateSmithWaterman();
                break;
            case 4:
                demonstrateNeedlemanWunsch();
                break;
            case 5:
                demonstrateFuzzySearch();
                break;
            case 6:
                demonstrateCompression();
                break;
            case 7:
                demonstrateClassicStringMatching();
                break;
            case 8:
                demonstrateAdvancedAlgorithms();
                break;
            case 9:
                demonstrateMCMCAndParallel();
                break;
            case 10:
                demonstrateWARPCTC();
                break;
            case 11:
                demonstrateAdvancedAlgorithms2();
                break;
            case 12:
                demonstrateHeuristicAndGraph();
                break;
            case 13:
                demonstrateVectorSimilarity();
                break;
            case 14:
                demonstrateVariantCalling();
                break;
            case 15:
                demonstrateMetagenomics();
                break;
            case 16:
                demonstrateChordHashing();
                break;
            case 17:
                demonstrateProbabilisticML();
                break;
            case 18:
                demonstrateConcurrentMultiSearch();
                break;
            case 19:
                demonstrateSkipGraph();
                break;
            case 20:
                demonstrateDancingLinks();
                break;
            case 21:
                demonstrateComparison();
                break;
            case 22:
                demonstrateExactMatch();
                demonstrateNaiveSearch();
                demonstrateSmithWaterman();
                demonstrateNeedlemanWunsch();
                demonstrateFuzzySearch();
                demonstrateCompression();
                demonstrateClassicStringMatching();
                demonstrateAdvancedAlgorithms();
                demonstrateMCMCAndParallel();
                demonstrateWARPCTC();
                demonstrateAdvancedAlgorithms2();
                demonstrateHeuristicAndGraph();
                demonstrateVectorSimilarity();
                demonstrateVariantCalling();
                demonstrateMetagenomics();
                demonstrateChordHashing();
                demonstrateProbabilisticML();
                demonstrateConcurrentMultiSearch();
                demonstrateSkipGraph();
                demonstrateDancingLinks();
                demonstrateComparison();
                break;
            case 0:
                cout << "Exiting..." << endl;
                break;
            default:
                cout << "Invalid option. Please try again." << endl;
        }
    } while (choice != 0);
    
    return 0;
}

