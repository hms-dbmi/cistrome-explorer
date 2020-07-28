import d3 from "./d3.js";
import range from "lodash/range";

export const HIGHLIGHTING_COLOR = ["#999", "gold", "#399AB6"][0];
export const HIGHLIGHTING_STROKE = ["black", "brown", "#265D8A"][0];
export const HIGHLIGHTING_OPACITY = 0.26;

export function drawRowHighlightRect(two, selectedRows, highlitRows, width, height) {
    if(!highlitRows || !selectedRows) {
        return;
    }
    const yScale = d3.scaleBand()
        .domain(range(selectedRows.length))
        .range([0, height]);
    const rowHeight = yScale.bandwidth();

    let aggregatedRows = 1;
    const sortedHighlitRows = highlitRows.slice().sort((a, b) => selectedRows.indexOf(a) - selectedRows.indexOf(b));
    sortedHighlitRows.forEach((d, i) => {
        const startY = yScale(selectedRows.indexOf(d));
        
        if(
            sortedHighlitRows.length > i + 1 && 
            Math.abs(selectedRows.indexOf(sortedHighlitRows[i + 1]) - selectedRows.indexOf(d)) === 1
        ) {
            // Aggregate highlighting rows for drawing the stroke only once
            aggregatedRows++;
            return;
        }

        const rect = two.makeRect(0, startY - rowHeight * (aggregatedRows - 1), width, rowHeight * aggregatedRows);
        rect.fill = HIGHLIGHTING_COLOR;
        rect.opacity = HIGHLIGHTING_OPACITY;
        rect.stroke = null; //HIGHLIGHTING_STROKE;

        aggregatedRows = 1;
    });
}