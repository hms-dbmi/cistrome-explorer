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
    const width = 140;
    const xMarginInitial = 0;
    const fontSize = 10;

    let left, textAlign;
    if(rowLinkPosition === "left") {
        left = x1 - xMarginInitial - width;
        textAlign = "right";
    } else if(rowLinkPosition === "right") {
        left = x0 + x1 + xMarginInitial;
        textAlign = "left";
    }

    // Scales
    const yScale = vega_scaleBand()
        .domain(range(rowInfo.length))
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
            {rowInfo.map((info, i) => (
                <div
                    key={i}
                    className="row-link"
                    style={{
                        height: `${rowHeight}px`,
                        width: `${width}px`
                    }}
                >
                    {info[rowLinkAttribute] ? (
                        <a 
                            href={info[rowLinkAttribute]}
                            title={info[rowLinkAttribute]}
                            target="_blank"
                            style={{
                                fontSize: `${fontSize}px`,
                                lineHeight: `${rowHeight}px`,
                                textAlign: textAlign,
                                [textAlign]: 0
                            }}
                        >
                            {info[rowLinkAttribute]}
                        </a>
                    ) : null}
                </div>
            ))}
        </div>
    );
}