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
        isLeft
    } = props;

    // Layouts and styles
    const isNominal = attribute.type === "nominal";
    const barAreaWidth = isNominal ? 20 : 50;
    const textAreaWidth = isNominal ? 50 : 20;
    const margin = 5;
    const titleFontSize = 12;
    const fontSize = 10;
    
    // Scales
    const valueExtent = [0, d3.extent(rowInfo.map(d => d[attribute.name]))[1]];   // Zero baseline
    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);
    const xScale = d3.scaleLinear()
        .domain(valueExtent)
        .range([0, barAreaWidth])
    const colorScale = attribute.type === "nominal" ? 
        d3.scaleOrdinal()
            .domain(Array.from(new Set(rowInfo.map(d => d[attribute.name]))))
            .range(d3.schemeSet3) : 
        d3.scaleLinear()
            .domain(valueExtent)
            .range([0, 1]);
    const rowHeight = yScale.bandwidth();

    // Render visual components for each row (i.e., bars and texts).
    const textAlign = isLeft ? "end" : "start";

    rowInfo.forEach((d, i) => {
        const barWidth = isNominal ? barAreaWidth : xScale(d[attribute.name]);
        
        const barLeft = left + (isLeft ? width - barWidth : 0);
        const textLeft = left + (isLeft ? width - barWidth - margin : barWidth + margin);
        const color = attribute.type === "nominal" ?
            colorScale(d[attribute.name]) : 
            d3.interpolateViridis(colorScale(d[attribute.name]));

        const rect = two.makeRect(barLeft, yScale(i), barWidth, rowHeight);
        rect.fill = color;

        // Render text labels when the space is enough.
        if(rowHeight >= fontSize){
            const text = two.makeText(textLeft, yScale(i) + rowHeight/2, textAreaWidth - 6, rowHeight, d[attribute.name])
            text.fill = d3.hsl(color).darker(3);
            text.fontsize = fontSize;
            text.align = textAlign;
            text.baseline = "middle";
            text.overflow = "ellipsis";

            console.log(textAreaWidth, d[attribute.name])
        }
    });

    // Title
    const titleLeft = left + (isLeft ? margin : width - margin);
    const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;

    // Draw a title of each dimension
    const titleText = attribute.name;
    const title = two.makeText(titleLeft, top, rowHeight, barAreaWidth, titleText);
    title.fill = "#9A9A9A";
    title.fontsize = titleFontSize;
    title.align = isLeft ? "end" : "start";
    title.baseline = "top";
    title.rotation = titleRotate;
}