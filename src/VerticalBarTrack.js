import range from 'lodash/range';
import d3 from './utils/d3.js';

/**
 * Component for rendering each attribute.
 * @prop {object} two An instance of the Two class.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} attribute The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 */
export function verticalBarTrack(props) {
    const {
        two,
        left, top, width, height,
        rowInfo, 
        attribute,
        isLeft,
        isCanvas
    } = props;

    // Layouts and styles
    const colWidth = 10;    // width of stacked bars
    const xMargin = 60;     // width of text area
    const xGap = 5;         // gap between bars and text
    const titleFontSize = 12;
    const fontSize = 10;
    
    // Scales
    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);
    const colorScale = attribute.type === "nominal" ? 
        d3.scaleOrdinal()
            .domain(Array.from(new Set(rowInfo.map(d => d[attribute.name]))))
            .range(d3.schemeSet3) : 
        d3.scaleLinear()
            .domain(d3.extent(rowInfo.map(d => d[attribute.name])))
            .range([0, 1]);
    const rowHeight = yScale.bandwidth();

    // Render

    // Title
    const titleLeft = left + (isLeft ? xGap : width - xGap);
    const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;

    // Draw a title of each dimension
    const titleText = `attribute: ${attribute.name} | type: ${attribute.type}`
    const title = two.makeText(titleLeft, top, rowHeight, colWidth, titleText);
    title.fill = "#9A9A9A";
    title.fontsize = titleFontSize;
    title.align = isLeft ? "end" : "start";
    title.baseline = "top";
    title.rotation = titleRotate;

    // Render visual components for each row (i.e., bars and texts)
    const barLeft = left + (isLeft ? xMargin : 0);
    const labelLeft = left + (isLeft ? xMargin - xGap : colWidth + xGap);
    const labelAlign = isLeft ? "end" : "start";

    rowInfo.forEach((d, i) => {
        const color = attribute.type === "nominal" ?
            colorScale(d[attribute.name]) : 
            d3.interpolateViridis(colorScale(d[attribute.name]));

        const rect = two.makeRect(barLeft, yScale(i), colWidth, rowHeight);
        rect.fill = color;

        // Render text labels when the space is enough
        if(rowHeight >= fontSize){
            const text = two.makeText(labelLeft, yScale(i) + rowHeight/2, colWidth, rowHeight, d[attribute.name])
            text.fill = d3.hsl(color).darker(3);
            text.fontsize = fontSize;
            text.align = labelAlign;
            text.baseline = "middle";
        }
    });
}