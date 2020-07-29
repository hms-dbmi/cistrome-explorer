import d3 from "./d3.js";

export const HIGHLIGHTING_COLOR = ["#333", "gold", "#399AB6"][0];
export const HIGHLIGHTING_STROKE = ["black", "brown", "#265D8A"][0];
export const HIGHLIGHTING_OPACITY = 0.24;

/**
 * Draw highlighting effect with a background rectangle.
 * @param {object} two An object of two.js class
 * @param {array} selectedRows The array of index for selected rows.
 * @param {array} highlitRows The array of index for highlighting rows.
 * @param {number} width Size of a track along x-axis.
 * @param {number} height Size of a track along y-axis.
 * @param {boolean} styles.isStroke Should draw stroke? 
 */
export function drawRowHighlightRect(two, selectedRows, highlitRows, width, height, styles) {
    if(!highlitRows || !selectedRows) {
        return;
    }

    const yScale = d3.scaleBand()
        .domain(selectedRows)
        .range([0, height]);
    const rowHeight = yScale.bandwidth();

    let aggregatedRows = 1;
    const sortedHighlitRows = highlitRows.slice().sort((a, b) => selectedRows.indexOf(a) - selectedRows.indexOf(b));
    sortedHighlitRows.forEach((d, i) => {
        const startY = yScale(d);
        
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
        rect.stroke = styles?.isStroke ? HIGHLIGHTING_STROKE : null;

        aggregatedRows = 1;
    });
}