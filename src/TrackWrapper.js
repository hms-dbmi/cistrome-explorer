import React from 'react';

import CistromeGroupLabels from './CistromeGroupLabels.js';


export default function TrackWrapper(props) {
    const { x0, options, track } = props;

    const x1 = track.position[0];
    const y = track.position[1];
    const height = track.dimensions[1];
    const rowNames = track.tilesetInfo.row_infos.map(d => d.split("\t"));

    console.log("TrackWrapper.render");
    return (
        <div className="cistrome-hgw-track-wrapper">
            {options.rowInfoPosition !== "hidden" ? 
                (<CistromeGroupLabels 
                    rowNames={rowNames}
                    x={x0 + x1}
                    y={y}
                    height={height}
                />) : null}
        </div>
    );
}