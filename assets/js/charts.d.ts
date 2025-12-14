export interface ChartData {
    labels: string[];
    datasets: {
        label: string;
        data: number[];
        backgroundColor?: string | string[];
        borderColor?: string | string[];
    }[];
}
export declare class ChartManager {
    private _charts;
    createBarChart(canvasId: string, data: ChartData, _options?: any): void;
    private createSVGChart;
    createLineChart(_canvasId: string, _data: ChartData, _options?: any): void;
    createPieChart(_canvasId: string, _data: ChartData, _options?: any): void;
}
export default ChartManager;
//# sourceMappingURL=charts.d.ts.map