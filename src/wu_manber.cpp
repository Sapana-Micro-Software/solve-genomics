#include "wu_manber.h"
#include <algorithm>
#include <cctype>
#include <climits>

WuManber::WuManber(int block_size) : block_size_(block_size) {
    if (block_size_ < 1) block_size_ = 2;
}

void WuManber::buildTables(const std::vector<std::string>& patterns) {
    patterns_ = patterns;
    shift_table_.clear();
    hash_table_.clear();
    prefix_table_.clear();
    
    if (patterns_.empty()) {
        return;
    }
    
    // Find minimum pattern length
    int min_len = INT_MAX;
    for (const auto& pattern : patterns_) {
        if (static_cast<int>(pattern.length()) < min_len) {
            min_len = static_cast<int>(pattern.length());
        }
    }
    
    if (min_len < block_size_) {
        block_size_ = min_len;
    }
    
    // Initialize shift table with default shift
    int default_shift = min_len - block_size_ + 1;
    
    // Build shift and hash tables
    for (size_t i = 0; i < patterns_.size(); ++i) {
        const std::string& pattern = patterns_[i];
        
        if (static_cast<int>(pattern.length()) < block_size_) {
            continue;
        }
        
        // Process each block in pattern
        for (int j = 0; j <= static_cast<int>(pattern.length()) - block_size_; ++j) {
            std::string block = pattern.substr(j, block_size_);
            
            // Convert to uppercase
            for (char& c : block) {
                c = std::toupper(c);
            }
            
            // Calculate shift
            int shift = static_cast<int>(pattern.length()) - j - block_size_;
            if (shift == 0) shift = 1;  // Minimum shift of 1
            
            // Update shift table
            if (shift_table_.find(block) == shift_table_.end() || 
                shift_table_[block] > shift) {
                shift_table_[block] = shift;
            }
            
            // Add to hash table
            hash_table_[block].push_back(static_cast<int>(i));
        }
        
        // Store prefix for verification
        if (pattern.length() >= static_cast<size_t>(block_size_)) {
            prefix_table_.push_back(pattern.substr(0, block_size_));
        } else {
            prefix_table_.push_back(pattern);
        }
    }
    
    // Set default shift for blocks not in patterns
    // This is handled during search
}

WuManber::MultiPatternResult WuManber::search(const std::string& text) {
    MultiPatternResult result;
    
    if (text.empty() || patterns_.empty()) {
        return result;
    }
    
    int min_pattern_len = INT_MAX;
    for (const auto& pattern : patterns_) {
        if (static_cast<int>(pattern.length()) < min_pattern_len) {
            min_pattern_len = static_cast<int>(pattern.length());
        }
    }
    
    if (min_pattern_len < block_size_) {
        return result;
    }
    
    size_t text_len = text.length();
    size_t pos = static_cast<size_t>(min_pattern_len) - block_size_;
    
    while (pos < text_len) {
        std::string block = getBlock(text, pos);
        
        if (shift_table_.find(block) != shift_table_.end()) {
            // Potential match - check hash table
            if (hash_table_.find(block) != hash_table_.end()) {
                for (int pattern_idx : hash_table_[block]) {
                    size_t match_pos = pos - (patterns_[pattern_idx].length() - block_size_);
                    
                    if (match_pos < text_len && 
                        verifyMatch(text, match_pos, pattern_idx)) {
                        result.matches[patterns_[pattern_idx]].push_back(match_pos);
                        result.total_matches++;
                    }
                }
            }
            
            // Shift by table value
            pos += shift_table_[block];
        } else {
            // Shift by default amount
            pos += min_pattern_len - block_size_ + 1;
        }
    }
    
    return result;
}

SearchResult WuManber::searchSingle(const std::string& text, const std::string& pattern) {
    SearchResult result;
    
    buildTables({pattern});
    MultiPatternResult multi_result = search(text);
    
    if (multi_result.matches.find(pattern) != multi_result.matches.end()) {
        result.positions.clear();
        for (size_t pos : multi_result.matches[pattern]) {
            result.positions.push_back(static_cast<int>(pos));
        }
        result.count = static_cast<int>(result.positions.size());
    }
    
    return result;
}

std::string WuManber::getBlock(const std::string& text, size_t pos) {
    if (pos + block_size_ > text.length()) {
        return "";
    }
    
    std::string block = text.substr(pos, block_size_);
    for (char& c : block) {
        c = std::toupper(c);
    }
    
    return block;
}

int WuManber::calculateShift(const std::string& block, int pattern_index) {
    const std::string& pattern = patterns_[pattern_index];
    
    // Find rightmost occurrence of block in pattern
    for (int i = static_cast<int>(pattern.length()) - block_size_; i >= 0; --i) {
        std::string pattern_block = pattern.substr(i, block_size_);
        for (char& c : pattern_block) {
            c = std::toupper(c);
        }
        
        if (pattern_block == block) {
            return static_cast<int>(pattern.length()) - i - block_size_;
        }
    }
    
    return static_cast<int>(pattern.length()) - block_size_ + 1;
}

bool WuManber::verifyMatch(const std::string& text, size_t pos, int pattern_index) {
    const std::string& pattern = patterns_[pattern_index];
    
    if (pos + pattern.length() > text.length()) {
        return false;
    }
    
    for (size_t i = 0; i < pattern.length(); ++i) {
        if (std::toupper(text[pos + i]) != std::toupper(pattern[i])) {
            return false;
        }
    }
    
    return true;
}

