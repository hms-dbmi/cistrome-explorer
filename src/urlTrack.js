import range from 'lodash/range';
import d3 from './utils/d3.js';

import { margin } from './visualizationTrack.js';

/**
 * Dendrogram for rendering an attribute of hierarchical structure.
 * @prop {object} two An instance of the Two class.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 */
export function urlTrack(props) {
    const {
        two,
        left, top, width, height,
        rowInfo, 
        fieldInfo,
        isLeft
    } = props;

    // Data, layouts and styles
    const { field, title } = fieldInfo;
    const fontSize = 10;
    const textAlign = isLeft ? "end" : "start";

    // Scales
    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);
    const rowHeight = yScale.bandwidth();

    if(rowHeight < fontSize) {
        // Bail out if there is not enough height per row to render links.
        return null;
    }

    rowInfo.forEach((info, i) => {
        const textTop = yScale(i);
        const textLeft = isLeft ? left + width - margin : left + margin;
        const titleField = title ? title : field;

        const text = two.makeText(textLeft, textTop + rowHeight/2, width, rowHeight, info[titleField]);
        text.fill = "#23527C";
        text.fontsize = fontSize;
        text.align = textAlign;
        text.baseline = "middle";
        text.overflow = "ellipsis";
    })
}