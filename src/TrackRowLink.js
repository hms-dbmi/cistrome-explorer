import React from 'react';
import range from 'lodash/range';

import { vega_scaleBand } from './utils-scales.js';

import './TrackRowLink.scss';

/**
 * Component for rendering links to row info metadata webpages.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {string} rowLinkAttribute The attribute used to obtain a URL from a row info JSON object.
 * @prop {string} rowLinkPosition The value of the `rowLinkPosition` option.
 */
export default function TrackRowLink(props) {

    const {
        trackX, trackY, 
        trackWidth, trackHeight, 
        rowInfo, 
        rowLinkPosition, rowLinkAttribute
    } = props;

    // Dimensions
    const top = trackY;
    const width = 140;
    const height = trackHeight;
    const fontSize = 10;

    let left, textAlign;
    if(rowLinkPosition === "left") {
        left = trackX - width;
        textAlign = "right";
    } else if(rowLinkPosition === "right") {
        left = trackWidth + trackX;
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
            className="cistrome-hgw-child"
            style={{
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