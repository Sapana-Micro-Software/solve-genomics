// Chart utilities for visualizing benchmark data

export interface ChartData {
  labels: string[];
  datasets: {
    label: string;
    data: number[];
    backgroundColor?: string | string[];
    borderColor?: string | string[];
  }[];
}

export class ChartManager {
  private _charts: Map<string, any> = new Map();

  createBarChart(canvasId: string, data: ChartData, _options?: any): void {
    const canvas = document.getElementById(canvasId) as HTMLCanvasElement;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) {
      // Fallback to SVG if canvas not available
      this.createSVGChart(canvasId, data);
      return;
    }

    // Chart.js implementation would go here
    // For now, we'll create a simple SVG-based visualization
    this.createSVGChart(canvasId, data);
  }

  private createSVGChart(containerId: string, data: ChartData): void {
    const container = document.getElementById(containerId.replace('Chart', ''));
    if (!container) return;

    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('width', '100%');
    svg.setAttribute('height', '400');
    svg.setAttribute('viewBox', '0 0 800 400');

    const maxValue = Math.max(...data.datasets[0].data);
    const barWidth = 700 / data.labels.length;
    const scale = 350 / maxValue;

    data.labels.forEach((label, index) => {
      const value = data.datasets[0].data[index];
      const barHeight = value * scale;
      const x = 50 + index * barWidth;
      const y = 350 - barHeight;

      // Bar
      const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      rect.setAttribute('x', x.toString());
      rect.setAttribute('y', y.toString());
      rect.setAttribute('width', (barWidth - 10).toString());
      rect.setAttribute('height', barHeight.toString());
      rect.setAttribute('fill', '#667eea');
      rect.setAttribute('rx', '4');
      svg.appendChild(rect);

      // Label
      const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      text.setAttribute('x', (x + barWidth / 2 - 5).toString());
      text.setAttribute('y', '380');
      text.setAttribute('text-anchor', 'middle');
      text.setAttribute('font-size', '12');
      text.setAttribute('fill', '#333');
      text.textContent = label;
      svg.appendChild(text);

      // Value
      const valueText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      valueText.setAttribute('x', (x + barWidth / 2 - 5).toString());
      valueText.setAttribute('y', (y - 5).toString());
      valueText.setAttribute('text-anchor', 'middle');
      valueText.setAttribute('font-size', '12');
      valueText.setAttribute('fill', '#333');
      valueText.setAttribute('font-weight', 'bold');
      valueText.textContent = value.toString();
      svg.appendChild(valueText);
    });

    container.appendChild(svg);
  }

  createLineChart(_canvasId: string, _data: ChartData, _options?: any): void {
    // Line chart implementation
  }

  createPieChart(_canvasId: string, _data: ChartData, _options?: any): void {
    // Pie chart implementation
  }
}

export default ChartManager;

