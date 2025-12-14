#include "warp_ctc.h"
#include <algorithm>
#include <cmath>
#include <cctype>

WarpCTC::WarpCTC(double blank_probability) : blank_probability_(blank_probability) {
    if (blank_probability_ < 0.0) blank_probability_ = 0.0;
    if (blank_probability_ > 1.0) blank_probability_ = 1.0;
}

WarpCTC::CTCResult WarpCTC::search(const std::string& sequence,
                                   const std::string& pattern,
                                   int max_gaps) {
    CTCResult result;
    
    if (sequence.empty() || pattern.empty() || pattern.length() > sequence.length() + max_gaps) {
        return result;
    }
    
    // Slide pattern over sequence
    size_t max_start = sequence.length() - pattern.length() + max_gaps;
    
    for (size_t start = 0; start <= max_start && start < sequence.length(); ++start) {
        // Extract window (may be longer than pattern due to gaps)
        size_t window_end = std::min(start + pattern.length() + max_gaps, sequence.length());
        std::string window = sequence.substr(start, window_end - start);
        
        // Calculate CTC alignment score
        double ctc_score = calculateCTCLoss(window, pattern);
        
        if (ctc_score > 0.1) {  // Threshold for valid alignment
            result.positions.push_back(start);
            result.alignment_scores.push_back(ctc_score);
        }
    }
    
    // Find best alignment
    if (!result.positions.empty()) {
        auto max_it = std::max_element(result.alignment_scores.begin(), 
                                       result.alignment_scores.end());
        size_t best_idx = std::distance(result.alignment_scores.begin(), max_it);
        size_t best_pos = result.positions[best_idx];
        
        size_t window_end = std::min(best_pos + pattern.length() + max_gaps, sequence.length());
        std::string best_window = sequence.substr(best_pos, window_end - best_pos);
        
        // Get aligned sequences
        auto aligned = align(best_window, pattern);
        result.aligned_pattern = aligned.first;
        result.aligned_sequence = aligned.second;
        result.total_score = *max_it;
        result.num_matches = static_cast<int>(result.positions.size());
    }
    
    return result;
}

double WarpCTC::calculateCTCLoss(const std::string& sequence, const std::string& pattern) {
    if (sequence.empty() || pattern.empty()) {
        return 0.0;
    }
    
    // Use forward algorithm to compute probability
    return forwardAlgorithm(sequence, pattern);
}

double WarpCTC::forwardAlgorithm(const std::string& sequence, const std::string& pattern) {
    std::string extended = extendWithBlanks(pattern);
    int T = sequence.length();
    int S = extended.length();
    
    if (T == 0 || S == 0) {
        return 0.0;
    }
    
    // Forward probabilities: alpha[t][s] = probability of being in state s at time t
    std::vector<std::vector<double>> alpha(T, std::vector<double>(S, 0.0));
    
    // Initialize: can start at first two states (blank or first pattern char)
    alpha[0][0] = blank_probability_;
    if (S > 1) {
        alpha[0][1] = (1.0 - blank_probability_) * emissionProbability(sequence[0], extended[1]);
    }
    
    // Forward pass
    for (int t = 1; t < T; ++t) {
        for (int s = 0; s < S; ++s) {
            double sum = 0.0;
            
            // Can come from same state (repeat)
            if (s < S) {
                sum += alpha[t-1][s] * transitionProbability(s, s, extended) * 
                       emissionProbability(sequence[t], extended[s]);
            }
            
            // Can come from previous state
            if (s > 0) {
                sum += alpha[t-1][s-1] * transitionProbability(s-1, s, extended) * 
                       emissionProbability(sequence[t], extended[s]);
            }
            
            // Can skip blank (if previous was blank and current is not)
            if (s > 1 && extended[s-1] == ' ') {
                sum += alpha[t-1][s-2] * transitionProbability(s-2, s, extended) * 
                       emissionProbability(sequence[t], extended[s]);
            }
            
            alpha[t][s] = sum;
        }
    }
    
    // Sum over all final states
    double total_prob = 0.0;
    for (int s = 0; s < S; ++s) {
        total_prob += alpha[T-1][s];
    }
    
    return total_prob;
}

double WarpCTC::backwardAlgorithm(const std::string& sequence, const std::string& pattern) {
    std::string extended = extendWithBlanks(pattern);
    int T = sequence.length();
    int S = extended.length();
    
    if (T == 0 || S == 0) {
        return 0.0;
    }
    
    // Backward probabilities: beta[t][s]
    std::vector<std::vector<double>> beta(T, std::vector<double>(S, 0.0));
    
    // Initialize: can end at any state
    for (int s = 0; s < S; ++s) {
        beta[T-1][s] = 1.0;
    }
    
    // Backward pass
    for (int t = T-2; t >= 0; --t) {
        for (int s = 0; s < S; ++s) {
            double sum = 0.0;
            
            // Can go to same state
            if (s < S) {
                sum += beta[t+1][s] * transitionProbability(s, s, extended) * 
                       emissionProbability(sequence[t+1], extended[s]);
            }
            
            // Can go to next state
            if (s < S-1) {
                sum += beta[t+1][s+1] * transitionProbability(s, s+1, extended) * 
                       emissionProbability(sequence[t+1], extended[s+1]);
            }
            
            // Can skip blank
            if (s < S-2 && extended[s+1] == ' ') {
                sum += beta[t+1][s+2] * transitionProbability(s, s+2, extended) * 
                       emissionProbability(sequence[t+1], extended[s+2]);
            }
            
            beta[t][s] = sum;
        }
    }
    
    // Sum over all initial states
    double total_prob = 0.0;
    for (int s = 0; s < std::min(2, S); ++s) {
        double init_prob = (s == 0) ? blank_probability_ : (1.0 - blank_probability_);
        total_prob += init_prob * beta[0][s] * emissionProbability(sequence[0], extended[s]);
    }
    
    return total_prob;
}

std::string WarpCTC::viterbiDecode(const std::string& sequence, const std::string& pattern) {
    std::string extended = extendWithBlanks(pattern);
    int T = sequence.length();
    int S = extended.length();
    
    if (T == 0 || S == 0) {
        return "";
    }
    
    // Viterbi: best path probability and backpointer
    std::vector<std::vector<double>> viterbi(T, std::vector<double>(S, 0.0));
    std::vector<std::vector<int>> backpointer(T, std::vector<int>(S, -1));
    
    // Initialize
    viterbi[0][0] = blank_probability_;
    if (S > 1) {
        viterbi[0][1] = (1.0 - blank_probability_) * emissionProbability(sequence[0], extended[1]);
    }
    
    // Viterbi recursion
    for (int t = 1; t < T; ++t) {
        for (int s = 0; s < S; ++s) {
            double best_prob = 0.0;
            int best_prev = -1;
            
            // From same state
            if (s < S && viterbi[t-1][s] > 0) {
                double prob = viterbi[t-1][s] * transitionProbability(s, s, extended) * 
                             emissionProbability(sequence[t], extended[s]);
                if (prob > best_prob) {
                    best_prob = prob;
                    best_prev = s;
                }
            }
            
            // From previous state
            if (s > 0 && viterbi[t-1][s-1] > 0) {
                double prob = viterbi[t-1][s-1] * transitionProbability(s-1, s, extended) * 
                             emissionProbability(sequence[t], extended[s]);
                if (prob > best_prob) {
                    best_prob = prob;
                    best_prev = s-1;
                }
            }
            
            // Skip blank
            if (s > 1 && extended[s-1] == ' ' && viterbi[t-1][s-2] > 0) {
                double prob = viterbi[t-1][s-2] * transitionProbability(s-2, s, extended) * 
                             emissionProbability(sequence[t], extended[s]);
                if (prob > best_prob) {
                    best_prob = prob;
                    best_prev = s-2;
                }
            }
            
            viterbi[t][s] = best_prob;
            backpointer[t][s] = best_prev;
        }
    }
    
    // Find best final state
    int best_final = 0;
    for (int s = 1; s < S; ++s) {
        if (viterbi[T-1][s] > viterbi[T-1][best_final]) {
            best_final = s;
        }
    }
    
    // Trace back
    std::string decoded;
    int state = best_final;
    for (int t = T-1; t >= 0; --t) {
        if (state >= 0 && state < S && extended[state] != ' ') {
            decoded = extended[state] + decoded;
        }
        if (t > 0) {
            state = backpointer[t][state];
        }
    }
    
    return decoded;
}

std::vector<std::pair<std::string, double>> WarpCTC::beamSearch(const std::string& sequence,
                                                                const std::string& pattern,
                                                                int beam_width) {
    std::vector<std::pair<std::string, double>> results;
    
    if (sequence.empty() || pattern.empty()) {
        return results;
    }
    
    std::string extended = extendWithBlanks(pattern);
    int T = sequence.length();
    int S = extended.length();
    
    // Beam: (state, probability, path)
    std::vector<std::tuple<int, double, std::string>> beam;
    beam.push_back({0, blank_probability_, ""});
    if (S > 1) {
        double prob = (1.0 - blank_probability_) * emissionProbability(sequence[0], extended[1]);
        beam.push_back({1, prob, std::string(1, extended[1])});
    }
    
    // Beam search
    for (int t = 1; t < T; ++t) {
        std::vector<std::tuple<int, double, std::string>> next_beam;
        
        for (const auto& [state, prob, path] : beam) {
            // Stay in same state
            if (state < S) {
                double new_prob = prob * transitionProbability(state, state, extended) * 
                                 emissionProbability(sequence[t], extended[state]);
                next_beam.push_back({state, new_prob, path});
            }
            
            // Move to next state
            if (state < S-1) {
                double new_prob = prob * transitionProbability(state, state+1, extended) * 
                                 emissionProbability(sequence[t], extended[state+1]);
                std::string new_path = path;
                if (extended[state+1] != ' ') {
                    new_path += extended[state+1];
                }
                next_beam.push_back({state+1, new_prob, new_path});
            }
            
            // Skip blank
            if (state < S-2 && extended[state+1] == ' ') {
                double new_prob = prob * transitionProbability(state, state+2, extended) * 
                                 emissionProbability(sequence[t], extended[state+2]);
                std::string new_path = path;
                if (extended[state+2] != ' ') {
                    new_path += extended[state+2];
                }
                next_beam.push_back({state+2, new_prob, new_path});
            }
        }
        
        // Keep top beam_width
        std::sort(next_beam.begin(), next_beam.end(),
            [](const std::tuple<int, double, std::string>& a,
               const std::tuple<int, double, std::string>& b) {
                return std::get<1>(a) > std::get<1>(b);
            });
        
        if (static_cast<int>(next_beam.size()) > beam_width) {
            next_beam.resize(beam_width);
        }
        
        beam = next_beam;
    }
    
    // Extract results
    for (const auto& [state, prob, path] : beam) {
        results.push_back({path, prob});
    }
    
    return results;
}

std::pair<std::string, std::string> WarpCTC::align(const std::string& sequence,
                                                   const std::string& pattern) {
    std::vector<std::vector<double>> matrix = buildAlignmentMatrix(sequence, pattern);
    std::string extended = extendWithBlanks(pattern);
    return tracebackAlignment(matrix, sequence, extended);
}

std::string WarpCTC::extendWithBlanks(const std::string& pattern) {
    std::string extended;
    extended.reserve(pattern.length() * 2 + 1);
    
    extended += ' ';  // Start with blank
    for (char c : pattern) {
        extended += std::toupper(c);
        extended += ' ';  // Blank between each character
    }
    
    return extended;
}

double WarpCTC::emissionProbability(char seq_char, char pattern_char) {
    seq_char = std::toupper(seq_char);
    pattern_char = std::toupper(pattern_char);
    
    if (pattern_char == ' ') {
        // Blank can emit any character with low probability
        return 0.1;
    }
    
    if (seq_char == pattern_char) {
        return 0.9;  // High probability for match
    }
    
    return 0.05;  // Low probability for mismatch
}

double WarpCTC::transitionProbability(int from_state, int to_state, const std::string& extended_pattern) {
    if (from_state < 0 || to_state < 0 || 
        from_state >= static_cast<int>(extended_pattern.length()) ||
        to_state >= static_cast<int>(extended_pattern.length())) {
        return 0.0;
    }
    
    // Can stay in same state (repeat)
    if (from_state == to_state) {
        return 0.5;
    }
    
    // Can move to next state
    if (to_state == from_state + 1) {
        return 0.4;
    }
    
    // Can skip blank (move two states if middle is blank)
    if (to_state == from_state + 2 && 
        from_state + 1 < static_cast<int>(extended_pattern.length()) &&
        extended_pattern[from_state + 1] == ' ') {
        return 0.1;
    }
    
    return 0.0;
}

std::vector<std::vector<double>> WarpCTC::buildAlignmentMatrix(const std::string& sequence,
                                                                const std::string& pattern) {
    std::string extended = extendWithBlanks(pattern);
    int T = sequence.length();
    int S = extended.length();
    
    std::vector<std::vector<double>> matrix(T, std::vector<double>(S, 0.0));
    
    // Initialize
    matrix[0][0] = blank_probability_;
    if (S > 1) {
        matrix[0][1] = (1.0 - blank_probability_) * emissionProbability(sequence[0], extended[1]);
    }
    
    // Fill matrix
    for (int t = 1; t < T; ++t) {
        for (int s = 0; s < S; ++s) {
            double sum = 0.0;
            
            if (s < S) {
                sum += matrix[t-1][s] * transitionProbability(s, s, extended) * 
                       emissionProbability(sequence[t], extended[s]);
            }
            
            if (s > 0) {
                sum += matrix[t-1][s-1] * transitionProbability(s-1, s, extended) * 
                       emissionProbability(sequence[t], extended[s]);
            }
            
            if (s > 1 && extended[s-1] == ' ') {
                sum += matrix[t-1][s-2] * transitionProbability(s-2, s, extended) * 
                       emissionProbability(sequence[t], extended[s]);
            }
            
            matrix[t][s] = sum;
        }
    }
    
    return matrix;
}

std::pair<std::string, std::string> WarpCTC::tracebackAlignment(
    const std::vector<std::vector<double>>& matrix,
    const std::string& sequence,
    const std::string& extended_pattern) {
    
    int T = matrix.size();
    int S = (T > 0) ? matrix[0].size() : 0;
    
    if (T == 0 || S == 0) {
        return {"", ""};
    }
    
    // Find best final state
    int best_state = 0;
    for (int s = 1; s < S; ++s) {
        if (matrix[T-1][s] > matrix[T-1][best_state]) {
            best_state = s;
        }
    }
    
    // Trace back
    std::string aligned_pattern;
    std::string aligned_sequence;
    int state = best_state;
    
    for (int t = T-1; t >= 0; --t) {
        if (state >= 0 && state < S) {
            if (extended_pattern[state] != ' ') {
                aligned_pattern = extended_pattern[state] + aligned_pattern;
                aligned_sequence = sequence[t] + aligned_sequence;
            } else {
                aligned_pattern = '-' + aligned_pattern;
                aligned_sequence = sequence[t] + aligned_sequence;
            }
        }
        
        // Find best previous state
        if (t > 0) {
            int best_prev = state;
            double best_prob = 0.0;
            
            if (state < S && matrix[t-1][state] > best_prob) {
                best_prob = matrix[t-1][state];
                best_prev = state;
            }
            
            if (state > 0 && matrix[t-1][state-1] > best_prob) {
                best_prob = matrix[t-1][state-1];
                best_prev = state - 1;
            }
            
            if (state > 1 && extended_pattern[state-1] == ' ' && 
                matrix[t-1][state-2] > best_prob) {
                best_prev = state - 2;
            }
            
            state = best_prev;
        }
    }
    
    return {aligned_pattern, aligned_sequence};
}

std::vector<std::pair<size_t, double>> WarpCTC::findAllAlignments(
    const std::string& sequence,
    const std::string& pattern) {
    
    std::vector<std::pair<size_t, double>> alignments;
    
    if (sequence.empty() || pattern.empty()) {
        return alignments;
    }
    
    // Slide pattern over sequence
    for (size_t start = 0; start <= sequence.length() - pattern.length(); ++start) {
        size_t end = std::min(start + pattern.length() * 2, sequence.length());
        std::string window = sequence.substr(start, end - start);
        
        double score = calculateCTCLoss(window, pattern);
        
        if (score > 0.1) {  // Threshold
            alignments.push_back({start, score});
        }
    }
    
    return alignments;
}

