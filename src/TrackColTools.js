import React from 'react';

import './TrackColTools.scss';

/**
 * Component for rendering genome interval selection tools.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {string} colToolsPosition The value of the `colToolsPosition` option.
 */
export default function TrackColTools(props) {

    const {
        trackX, trackY, 
        trackWidth, trackHeight,
        colToolsPosition,
        onSelectGenomicInterval
    } = props;

    const left = trackX;
    const width = trackWidth;
    const height = 30;

    let top;
    if(colToolsPosition === "top") {
        top = trackY - height;
    } else if(colToolsPosition === "bottom") {
        top = trackY + trackHeight;
    }

    return (
        <div
            style={{
                position: 'absolute',
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height}px`
            }}
        >
            <div className="col-tools">
                <button onClick={onSelectGenomicInterval}>Select genomic interval</button>
            </div>
        </div>
    );
};