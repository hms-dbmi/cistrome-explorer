import range from 'lodash/range';
import d3 from './utils/d3.js';

import { margin } from './visualizationTrack.js';

/**
 * Bar chart for rendering an attribute.
 * @prop {object} two An instance of the Two class.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 */
export function verticalBarTrack(props) {
    const {
        two,
        left, top, width, height,
        rowInfo, 
        fieldInfo,
        isLeft
    } = props;

    // Data, layouts and styles
    const { field, type } = fieldInfo;
    const isNominal = type === "nominal";
    const barAreaWidth = isNominal ? 20 : width - 20;
    const textAreaWidth = isNominal ? 50 : 20;
    const fontSize = 10;
    
    // Scales
    const valueExtent = [0, d3.extent(rowInfo.map(d => d[field]))[1]];   // Zero baseline
    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);
    const xScale = d3.scaleLinear()
        .domain(valueExtent)
        .range([0, barAreaWidth])
    const colorScale = isNominal ? 
        d3.scaleOrdinal()
            .domain(Array.from(new Set(rowInfo.map(d => d[field]))))
            .range(d3.schemeSet3) : 
        d3.scaleLinear()
            .domain(valueExtent)
            .range([0, 1]);
    const rowHeight = yScale.bandwidth();

    // Render visual components for each row (i.e., bars and texts).
    const textAlign = isLeft ? "end" : "start";
    let aggregateStartIdx = -1, sameCategoriesNearby = 1;

    rowInfo.forEach((d, i) => {

        // To aggregate bars, check if there is a same category on the next row.
        if(i + 1 < rowInfo.length && d[field] === rowInfo[i+1][field]) {
            if(aggregateStartIdx === -1) {
                aggregateStartIdx = i;
            }
            sameCategoriesNearby++;
            return;
        }

        const barTop = aggregateStartIdx !== -1 ? yScale(aggregateStartIdx) : yScale(i);
        const barHeight = rowHeight * sameCategoriesNearby;
        const barWidth = isNominal ? barAreaWidth : xScale(d[field]);        
        const barLeft = left + (isLeft ? width - barWidth : 0);
        const textLeft = left + (isLeft ? width - barWidth - margin : barWidth + margin);
        const color = isNominal ? colorScale(d[field]) : 
            d3.interpolateViridis(colorScale(d[field]));

        const rect = two.makeRect(barLeft, barTop, barWidth, barHeight);
        rect.fill = color;

        // Render text labels when the space is enough.
        if(barHeight >= fontSize){
            const text = two.makeText(textLeft, barTop + barHeight/2, textAreaWidth, barHeight, d[field]);
            text.fill = d3.hsl(color).darker(3);
            text.fontsize = fontSize;
            text.align = textAlign;
            text.baseline = "middle";
            text.overflow = "ellipsis";
        }

        aggregateStartIdx = -1;
        sameCategoriesNearby = 1;
    });
}