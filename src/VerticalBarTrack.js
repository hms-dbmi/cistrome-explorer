import range from 'lodash/range';
import d3 from './d3.js';

/**
 * Component for rendering each attribute.
 * @prop {object} ref The context of canvas or svg to render this view.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} attribute The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {boolean} isCanvas Is this view rendered using Canvas?
 */
export function verticalBarTrack(props) {
    const {
        ref,
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
    const fontFamily = "Arial";
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
    if(isCanvas) {
        // Title
        const titleLeft = left + (isLeft ? xGap : width - xGap);
        const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;

        ref.fillStyle = "#9A9A9A";
        ref.font = `${titleFontSize}px ${fontFamily}`;
        ref.textAlign = isLeft ? "end" : "start";
        ref.textBaseline = "top";
        ref.translate(titleLeft, top);
        ref.rotate(titleRotate);
        ref.fillText(`attribute: ${attribute.name} | type: ${attribute.type}`, 0, top);
        ref.rotate(-titleRotate);
        ref.translate(-titleLeft, top);

        // Render visual components for each row (i.e., bars and texts)
        const barLeft = left + (isLeft ? xMargin : 0);
        const labelLeft = left + (isLeft ? xMargin - xGap : colWidth + xGap);
        const labelAlign = isLeft ? "end" : "start";

        rowInfo.forEach((d, i) => {
            const color = attribute.type === "nominal" ?
             colorScale(d[attribute.name]) : 
             d3.interpolateViridis(colorScale(d[attribute.name]));

            ref.fillStyle = color;
            ref.fillRect(barLeft, yScale(i), colWidth, rowHeight);

            // Render text labels when the space is enough
            if(rowHeight >= fontSize){
                ref.fillStyle = d3.hsl(color).darker(3);
                ref.font = `${fontSize}px ${fontFamily}`;
                ref.textAlign = labelAlign;
                ref.textBaseline = "middle";
                ref.fillText(d[attribute.name], labelLeft, yScale(i) + rowHeight / 2.0);
            }
        });
    }
}