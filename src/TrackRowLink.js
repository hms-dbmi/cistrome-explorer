import React from 'react';
import range from 'lodash/range';

import { vega_scaleBand } from './utils-scales.js';

import './TrackRowLink.scss';

/**
 * Component for rendering links to row info metadata webpages.
 * @prop {number} x0
 * @prop {number} x1
 * @prop {number} y0
 * @prop {number} y1
 * @prop {number} height The track height.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {string} rowLinkAttribute The attribute used to obtain a URL from a row info JSON object.
 * @prop {string} rowLinkPosition The value of the `rowLinkPosition` option.
 */
export default function TrackRowLink(props) {

    const {
        x0, x1, y1, height, rowInfo, 
        rowLinkPosition, rowLinkAttribute
    } = props;

    // Dimensions
    const top = y1;
    const width = 120;
    const xMarginInitial = 5;
    const fontSize = 9;

    let left;
    if(rowLinkPosition === "left") {
        left = x1 - xMarginInitial - width;
    } else if(rowLinkPosition === "right") {
        left = x0 + x1 + xMarginInitial;
    }

    // Scales
    const yScale = vega_scaleBand()
        .domain(range(rowInfo.length - 1))
        .range([0, height]);

    const rowHeight = yScale.bandwidth();
    if(rowHeight < fontSize) {
        // Bail out if there is not enough height per row to render links.
        return null;
    }

    return (
        <div
            style={{
                position: 'absolute',
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height}px`,
            }}
        >
            <p>TODO</p>
        </div>
    );
}