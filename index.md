---
layout: default
title: Solve Genomics
description: Comprehensive DNA Sequence Alignment and Pattern Matching Algorithms
---

<div class="gradient-header relative overflow-hidden">
    <div class="absolute inset-0 opacity-10">
        <img src="{{ '/assets/images/svg/dna-helix.svg' | relative_url }}" alt="" class="w-full h-full object-cover">
    </div>
    <div class="container mx-auto px-4 relative z-10">
        <div class="flex flex-col md:flex-row items-center justify-between">
            <div class="flex-1 text-center md:text-left mb-8 md:mb-0">
                <h1 class="text-5xl font-bold mb-4">Solve Genomics</h1>
                <p class="text-xl opacity-95 max-w-3xl mx-auto md:mx-0">
                    Comprehensive DNA Sequence Alignment and Pattern Matching Algorithms
                </p>
                <div class="mt-8 text-lg opacity-90">
                    <strong>Shyamal Suhana Chandra</strong><br>
                    <span class="text-lg">Sapana Micro Software</span><br>
                    <span class="text-base">Implementation, Analysis, and Performance Evaluation</span>
                </div>
            </div>
            <div class="flex-1 flex justify-center md:justify-end">
                <img src="{{ '/assets/images/svg/dna-helix.svg' | relative_url }}" alt="DNA Helix" class="w-64 h-64 svg-glow svg-animate">
            </div>
        </div>
    </div>
</div>

<section id="overview" class="py-20">
    <div class="container mx-auto px-4">
        <h2 class="section-title">
            We present a comprehensive implementation of 25+ DNA sequence alignment and pattern matching algorithms with detailed benchmarks and analysis.
        </h2>
        
        <div class="card max-w-4xl mx-auto mb-12">
            <h2 class="text-3xl font-bold text-primary-600 mb-6">Abstract</h2>
            <p class="text-lg leading-relaxed">
                This project presents a comprehensive implementation and analysis of DNA sequence alignment and pattern matching algorithms. 
                We implement and compare multiple approaches including exact matching, approximate matching, dynamic programming algorithms 
                (Smith-Waterman, Needleman-Wunsch), fuzzy search with edit distance, classic string matching algorithms (Rabin-Karp, KMP, 
                Boyer-Moore), compression techniques (grammar-based and lossy), embedding-based search, deep learning approaches (CNN, 
                lightweight transformers), MCMC-based pattern evolution, WARP-CTC alignment, parallel/distributed search methods, 
                concurrent multi-technique search, skip-graph hierarchical indexing, and dancing links for exact cover problems. 
                We provide detailed performance benchmarks on sequences with varying complexity (high entropy vs. high repetition) and 
                analyze scalability characteristics. Our results demonstrate the trade-offs between accuracy, speed, and memory usage 
                across different algorithm classes, providing guidance for algorithm selection based on use case requirements.
            </p>
        </div>

        <div class="grid grid-cols-1 md:grid-cols-4 gap-6 mb-12">
            <div class="card text-center relative overflow-hidden group">
                <div class="absolute top-0 right-0 w-20 h-20 opacity-10">
                    <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="text-4xl font-bold text-primary-600 mb-2 relative z-10">25+</div>
                <div class="text-gray-600 relative z-10">Algorithms</div>
            </div>
            <div class="card text-center relative overflow-hidden group">
                <div class="absolute top-0 right-0 w-20 h-20 opacity-10">
                    <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="text-4xl font-bold text-primary-600 mb-2 relative z-10">7,000+</div>
                <div class="text-gray-600 relative z-10">Test Cases</div>
            </div>
            <div class="card text-center relative overflow-hidden group">
                <div class="absolute top-0 right-0 w-20 h-20 opacity-10">
                    <img src="{{ '/assets/images/svg/performance-chart.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="text-4xl font-bold text-primary-600 mb-2 relative z-10">100%</div>
                <div class="text-gray-600 relative z-10">Code Coverage</div>
            </div>
            <div class="card text-center relative overflow-hidden group">
                <div class="absolute top-0 right-0 w-20 h-20 opacity-10">
                    <img src="{{ '/assets/images/svg/dna-sequence.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="text-4xl font-bold text-primary-600 mb-2 relative z-10">C++</div>
                <div class="text-gray-600 relative z-10">Implementation</div>
            </div>
        </div>
    </div>
</section>

<section id="algorithms" class="py-20 bg-white">
    <div class="container mx-auto px-4">
        <h2 class="section-title">Algorithm Categories</h2>
        
        <div class="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6" id="algorithms-container">
            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/dna-sequence.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/dna-sequence.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Exact Matching</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Exact Match</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Naive Search</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Rabin-Karp</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> KMP</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Boyer-Moore</li>
                    </ul>
                </div>
            </div>

            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Approximate Matching</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Fuzzy Search</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Edit Distance</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Levenshtein</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Hamming</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> DNA-specific</li>
                    </ul>
                </div>
            </div>

            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/performance-chart.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/performance-chart.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Dynamic Programming</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Smith-Waterman</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Needleman-Wunsch</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Space-optimized</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Affine gap</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Banded alignment</li>
                    </ul>
                </div>
            </div>

            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Compression</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Grammar-based</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Lossy compression</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Association lists</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Pattern approximation</li>
                    </ul>
                </div>
            </div>

            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/dna-helix.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/dna-helix.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Modern ML</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Embedding search</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> CNN</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Lightweight LLM</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> LSTM/GRU</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Attention models</li>
                    </ul>
                </div>
            </div>

            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Advanced Methods</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> MCMC Evolution</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> DDMCMC</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> WARP-CTC</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Aho-Corasick</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Suffix Trees</li>
                    </ul>
                </div>
            </div>

            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Parallel & Distributed</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Multi-threaded</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Map-Reduce</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Work-stealing</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Pipeline</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Concurrent</li>
                    </ul>
                </div>
            </div>

            <div class="card relative overflow-hidden group hover:scale-105 transition-transform duration-300">
                <div class="absolute top-0 right-0 w-32 h-32 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/performance-chart.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="relative z-10">
                    <div class="flex items-center mb-4">
                        <img src="{{ '/assets/images/svg/performance-chart.svg' | relative_url }}" alt="" class="w-12 h-12 mr-3 svg-glow">
                        <h3 class="text-xl font-bold text-primary-600">Indexing Structures</h3>
                    </div>
                    <ul class="space-y-2 text-gray-700">
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Skip-Graph</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Chord DHT</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Dancing Links</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Wu-Manber</li>
                        <li class="flex items-center"><span class="text-primary-500 mr-2">✓</span> Bitap</li>
                    </ul>
                </div>
            </div>
        </div>
    </div>
</section>

<section id="downloads" class="py-20">
    <div class="container mx-auto px-4">
        <h2 class="section-title">Downloads</h2>
        
        <div class="card max-w-4xl mx-auto">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">Papers and Documentation</h2>
            
            <div class="space-y-4">
                <div class="flex justify-between items-center p-4 bg-gradient-to-r from-gray-50 to-gray-100 rounded-lg hover:shadow-lg transition-all duration-300 group relative overflow-hidden">
                    <div class="absolute top-0 right-0 w-24 h-24 opacity-5 group-hover:opacity-10 transition-opacity">
                        <img src="{{ '/assets/images/svg/dna-sequence.svg' | relative_url }}" alt="" class="w-full h-full">
                    </div>
                    <div class="flex items-center flex-1 relative z-10">
                        <img src="{{ '/assets/images/svg/performance-chart.svg' | relative_url }}" alt="" class="w-12 h-12 mr-4 svg-glow">
                        <div class="flex-1">
                            <a href="{{ '/docs/pdfs/paper.pdf' | relative_url }}" 
                               class="text-lg font-semibold text-primary-600 hover:underline flex items-center" download>
                                <span class="mr-2">📄</span> Paper (PDF)
                            </a>
                            <p class="text-sm text-gray-600 mt-1">Comprehensive paper with all algorithms, implementations, and analysis</p>
                        </div>
                    </div>
                    <span class="text-sm text-gray-500 relative z-10">~9 KB</span>
                </div>

                <div class="flex justify-between items-center p-4 bg-gradient-to-r from-gray-50 to-gray-100 rounded-lg hover:shadow-lg transition-all duration-300 group relative overflow-hidden">
                    <div class="absolute top-0 right-0 w-24 h-24 opacity-5 group-hover:opacity-10 transition-opacity">
                        <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="" class="w-full h-full">
                    </div>
                    <div class="flex items-center flex-1 relative z-10">
                        <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="" class="w-12 h-12 mr-4 svg-glow">
                        <div class="flex-1">
                            <a href="{{ '/docs/pdfs/presentation.pdf' | relative_url }}" 
                               class="text-lg font-semibold text-primary-600 hover:underline flex items-center" download>
                                <span class="mr-2">📊</span> Presentation (PDF)
                            </a>
                            <p class="text-sm text-gray-600 mt-1">Beamer presentation with visualizations and results</p>
                        </div>
                    </div>
                    <span class="text-sm text-gray-500 relative z-10">~235 KB</span>
                </div>

                <div class="flex justify-between items-center p-4 bg-gradient-to-r from-gray-50 to-gray-100 rounded-lg hover:shadow-lg transition-all duration-300 group relative overflow-hidden">
                    <div class="absolute top-0 right-0 w-24 h-24 opacity-5 group-hover:opacity-10 transition-opacity">
                        <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="" class="w-full h-full">
                    </div>
                    <div class="flex items-center flex-1 relative z-10">
                        <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="" class="w-12 h-12 mr-4 svg-glow">
                        <div class="flex-1">
                            <a href="{{ '/docs/pdfs/benchmark_results.pdf' | relative_url }}" 
                               class="text-lg font-semibold text-primary-600 hover:underline flex items-center" download>
                                <span class="mr-2">📈</span> Benchmark Results (PDF)
                            </a>
                            <p class="text-sm text-gray-600 mt-1">Comprehensive benchmark results with charts, graphs, and performance analysis</p>
                        </div>
                    </div>
                    <span class="text-sm text-gray-500 relative z-10">~172 KB</span>
                </div>
            </div>

            <h2 class="text-2xl font-bold text-primary-600 mb-6 mt-10">Source Code</h2>
            
            <div class="flex justify-between items-center p-4 bg-gradient-to-r from-gray-50 to-gray-100 rounded-lg hover:shadow-lg transition-all duration-300 group relative overflow-hidden">
                <div class="absolute top-0 right-0 w-24 h-24 opacity-5 group-hover:opacity-10 transition-opacity">
                    <img src="{{ '/assets/images/svg/dna-helix.svg' | relative_url }}" alt="" class="w-full h-full">
                </div>
                <div class="flex items-center flex-1 relative z-10">
                    <img src="{{ '/assets/images/svg/dna-helix.svg' | relative_url }}" alt="" class="w-12 h-12 mr-4 svg-glow">
                    <div class="flex-1">
                        <a href="https://github.com/Sapana-Micro-Software/solve-genomics" 
                           class="text-lg font-semibold text-primary-600 hover:underline flex items-center" target="_blank">
                            <span class="mr-2">🔗</span> GitHub Repository
                        </a>
                        <p class="text-sm text-gray-600 mt-1">Complete C++ implementation with tests and benchmarks</p>
                    </div>
                </div>
                <span class="text-sm text-gray-500 relative z-10">GitHub</span>
            </div>
        </div>
    </div>
</section>

<section id="code" class="py-20 bg-white">
    <div class="container mx-auto px-4">
        <h2 class="section-title">Code</h2>
        
        <div class="card max-w-4xl mx-auto">
            <h2 class="text-2xl font-bold text-primary-600 mb-4">Quick Start</h2>
            <p class="mb-4">Clone the repository and build the project:</p>
            <pre class="bg-gray-900 text-gray-100 p-4 rounded-lg overflow-x-auto"><code>git clone https://github.com/Sapana-Micro-Software/solve-genomics.git
cd solve-genomics
mkdir build && cd build
cmake ..
make -j4
./bin/dna_aligner</code></pre>

            <h2 class="text-2xl font-bold text-primary-600 mb-4 mt-8">Example Usage</h2>
            <pre class="bg-gray-900 text-gray-100 p-4 rounded-lg overflow-x-auto"><code># Exact match search
./bin/dna_aligner exact

# Fuzzy search with edit distance
./bin/dna_aligner fuzzy

# Concurrent multi-technique search
./bin/dna_aligner concurrent

# Skip-graph indexing
./bin/dna_aligner skip

# Dancing links exact cover
./bin/dna_aligner dancing</code></pre>

            <div class="mt-6">
                <a href="https://github.com/Sapana-Micro-Software/solve-genomics" 
                   class="btn mr-4" target="_blank">View on GitHub</a>
                <a href="https://github.com/Sapana-Micro-Software/solve-genomics/blob/main/README.md" 
                   class="btn btn-secondary" target="_blank">View README</a>
            </div>
        </div>
    </div>
</section>

<section id="results" class="py-20">
    <div class="container mx-auto px-4">
        <h2 class="section-title">Key Results</h2>
        
        <div class="card max-w-4xl mx-auto mb-8">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">Performance Highlights</h2>
            <ul class="space-y-3 text-lg">
                <li><strong>Fastest Exact Matching:</strong> Boyer-Moore (25μs) and KMP (38μs) for 1000-base sequences</li>
                <li><strong>Most Accurate:</strong> Dynamic programming algorithms (Smith-Waterman, Needleman-Wunsch) with optimal alignments</li>
                <li><strong>Best Scalability:</strong> Parallel work-stealing achieves 7.8x speedup with 8 threads</li>
                <li><strong>Best Compression:</strong> Grammar compression achieves 0.3x ratio for low-entropy sequences</li>
                <li><strong>Most Versatile:</strong> Concurrent multi-technique search provides comprehensive pattern matching</li>
                <li><strong>Best for Long Sequences:</strong> Skip-graph indexing provides O(1) hash table lookup</li>
                <li><strong>Best for Approximate:</strong> Fuzzy search with edit distance (95μs average)</li>
            </ul>
        </div>

        <div class="card max-w-4xl mx-auto">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">Algorithm Selection Guidelines</h2>
            <ul class="space-y-3 text-lg">
                <li>Use <strong>Boyer-Moore</strong> or <strong>KMP</strong> for exact pattern matching</li>
                <li>Use <strong>Smith-Waterman</strong> for local alignment with optimal accuracy</li>
                <li>Use <strong>Skip-Graph</strong> for indexed search on long sequences</li>
                <li>Use <strong>Concurrent Multi-Technique</strong> for comprehensive pattern matching</li>
                <li>Use <strong>Parallel Work-Stealing</strong> for large-scale distributed search</li>
                <li>Use <strong>Grammar Compression</strong> for storage of repetitive sequences</li>
                <li>Use <strong>Embedding Search</strong> for similarity-based retrieval</li>
            </ul>
        </div>

        <!-- Performance Chart -->
        <div class="card max-w-4xl mx-auto mt-8">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">Performance Comparison</h2>
            <div class="flex flex-col md:flex-row items-center gap-8">
                <div class="flex-1">
                    <img src="{{ '/assets/images/svg/performance-chart.svg' | relative_url }}" alt="Performance Chart" class="w-full h-auto svg-glow">
                </div>
                <div id="performanceChartContainer" class="flex-1 w-full h-96"></div>
            </div>
        </div>
        
        <!-- Algorithm Flow Visualization -->
        <div class="card max-w-4xl mx-auto mt-8">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">Algorithm Processing Flow</h2>
            <img src="{{ '/assets/images/svg/algorithm-flow.svg' | relative_url }}" alt="Algorithm Flow" class="w-full h-auto svg-glow svg-animate">
        </div>
        
        <!-- DNA Sequence Visualization -->
        <div class="card max-w-4xl mx-auto mt-8">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">DNA Sequence Pattern Matching</h2>
            <img src="{{ '/assets/images/svg/dna-sequence.svg' | relative_url }}" alt="DNA Sequence" class="w-full h-auto svg-glow">
        </div>
        
        <!-- Network Graph -->
        <div class="card max-w-4xl mx-auto mt-8">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">Algorithm Network Graph</h2>
            <img src="{{ '/assets/images/svg/network-graph.svg' | relative_url }}" alt="Network Graph" class="w-full h-auto svg-glow">
        </div>
    </div>
</section>

<section id="related-work" class="py-20 bg-white">
    <div class="container mx-auto px-4">
        <h2 class="section-title">Related Work</h2>
        
        <div class="card max-w-4xl mx-auto">
            <h2 class="text-2xl font-bold text-primary-600 mb-6">Recent Advances in Sequence Analysis</h2>
            
            <div class="mb-6">
                <h3 class="text-xl font-semibold mb-3">Protein Language Models (PLMs)</h3>
                <p class="text-gray-700 mb-3">
                    Protein language models have emerged as transformative tools for understanding and interpreting 
                    protein sequences, enabling advances in structure prediction, functional annotation, and variant 
                    effect assessment directly from sequence alone.
                </p>
                <p class="text-gray-700 mb-3">
                    Recent developments include <strong>Bag-of-Mer (BoM) pooling</strong>, a biologically inspired 
                    strategy for aggregating amino acid embeddings that captures both local motifs and long-range 
                    interactions, and <strong>ARIES</strong>, a highly scalable multiple-sequence alignment algorithm 
                    that leverages PLM embeddings to achieve superior accuracy even in low-identity regions.
                </p>
                <p class="text-sm text-gray-600 mt-4">
                    <strong>Reference:</strong> M. Singh, "Advancing protein sequence analysis with protein language models," 
                    <em>MIT CSAIL Bioinformatics Seminar</em>, December 10, 2025. 
                    <a href="https://www.csail.mit.edu/event/advancing-protein-sequence-analysis-protein-language-models" 
                       class="text-primary-600 hover:underline" target="_blank">View presentation</a>
                </p>
            </div>
        </div>
    </div>
</section>

<script src="{{ '/assets/js/main.js' | relative_url }}"></script>
<script>
    // Initialize charts when page loads
    document.addEventListener('DOMContentLoaded', () => {
        if (typeof ChartManager !== 'undefined') {
            const chartManager = new ChartManager();
            const chartData = {
                labels: ['Exact Match', 'Naive', 'Rabin-Karp', 'KMP', 'Boyer-Moore'],
                datasets: [{
                    label: 'Execution Time (μs)',
                    data: [45, 48, 52, 38, 25]
                }]
            };
            chartManager.createBarChart('performanceChartContainer', chartData);
        }
    });
</script>

