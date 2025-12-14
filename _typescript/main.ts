// Main TypeScript file for Solve Genomics website

interface Algorithm {
  name: string;
  category: string;
  complexity: string;
  description: string;
}

interface BenchmarkResult {
  algorithm: string;
  time: number;
  memory: number;
  accuracy: number;
}

class SolveGenomicsApp {
  private algorithms: Algorithm[] = [];
  private _benchmarks: BenchmarkResult[] = [];

  constructor() {
    this.init();
  }

  private init(): void {
    this.loadAlgorithms();
    this.setupEventListeners();
    this.initializeCharts();
    this.animateOnScroll();
  }

  private loadAlgorithms(): void {
    // Load algorithm data
    this.algorithms = [
      { name: "Boyer-Moore", category: "Exact Matching", complexity: "O(n/m)", description: "Fastest exact matching algorithm" },
      { name: "KMP", category: "Exact Matching", complexity: "O(n+m)", description: "Knuth-Morris-Pratt with failure function" },
      { name: "Smith-Waterman", category: "Dynamic Programming", complexity: "O(n*m)", description: "Local sequence alignment" },
      { name: "Skip-Graph", category: "Indexing", complexity: "O(log n)", description: "Hierarchical indexing structure" },
      { name: "Concurrent Multi-Technique", category: "Parallel", complexity: "O(n/p)", description: "Multiple algorithms in parallel" },
    ];
  }

  private setupEventListeners(): void {
    // Navigation smooth scroll
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
      anchor.addEventListener('click', (e: Event) => {
        e.preventDefault();
        const target = document.querySelector((anchor as HTMLAnchorElement).getAttribute('href') || '');
        if (target) {
          target.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }
      });
    });

    // Download buttons
    document.querySelectorAll('.download-btn').forEach(btn => {
      btn.addEventListener('click', (_e) => {
        const url = (btn as HTMLElement).dataset.url;
        if (url) {
          window.open(url, '_blank');
        }
      });
    });
  }

  private initializeCharts(): void {
    // Initialize performance charts using Chart.js
    const canvas = document.getElementById('performanceChart') as HTMLCanvasElement;
    if (canvas) {
      const ctx = canvas.getContext('2d');
      if (ctx) {
        // Chart will be initialized here
        this.createPerformanceChart(ctx);
      }
    }
  }

  private createPerformanceChart(_ctx: CanvasRenderingContext2D): void {
    // Performance data
    // Chart.js would be used here if included
    // const data = {
    //   labels: ['Exact Match', 'Naive', 'Rabin-Karp', 'KMP', 'Boyer-Moore'],
    //   datasets: [{
    //     label: 'Execution Time (μs)',
    //     data: [45, 48, 52, 38, 25],
    //     backgroundColor: 'rgba(102, 126, 234, 0.5)',
    //     borderColor: 'rgba(102, 126, 234, 1)',
    //     borderWidth: 2
    //   }]
    // };
    // new Chart(ctx, { type: 'bar', data: data, options: {...} });
  }

  private animateOnScroll(): void {
    const observer = new IntersectionObserver((entries) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          entry.target.classList.add('animate-fade-in');
        }
      });
    }, { threshold: 0.1 });

    document.querySelectorAll('.card, .feature-card').forEach(el => {
      observer.observe(el);
    });
  }

  public filterAlgorithms(category: string): void {
    const filtered = this.algorithms.filter(alg => alg.category === category);
    this.renderAlgorithms(filtered);
  }

  private renderAlgorithms(algorithms: Algorithm[]): void {
    const container = document.getElementById('algorithms-container');
    if (container) {
      container.innerHTML = algorithms.map(alg => `
        <div class="card">
          <h3 class="text-xl font-bold text-primary-600 mb-2">${alg.name}</h3>
          <p class="text-gray-600 mb-2">${alg.description}</p>
          <span class="text-sm text-gray-500">Complexity: ${alg.complexity}</span>
        </div>
      `).join('');
    }
  }
}

// Initialize app when DOM is ready
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', () => {
    new SolveGenomicsApp();
  });
} else {
  new SolveGenomicsApp();
}

// Export for module usage
export default SolveGenomicsApp;

