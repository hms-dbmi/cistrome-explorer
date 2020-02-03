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
        colToolsPosition
    } = props;

    const left = trackX;
    const top = trackY + trackHeight;
    const width = trackWidth;
    const height = 100;

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
                TODO: implement genome interval selection tools
            </div>
        </div>
    );
};