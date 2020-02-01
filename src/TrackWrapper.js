import React from 'react';

import TrackRowInfo from './TrackRowInfo.js';


export default function TrackWrapper(props) {
    const { x0, y0, options, track } = props;

    const x1 = track.position[0];
    const y1 = track.position[1];
    const height = track.dimensions[1];
    let rowInfo = [];
    try {
        rowInfo = track.tilesetInfo.row_infos.map(JSON.parse);
    } catch(e) {
        // Tileset info probably has not yet loaded:
        console.log(e);
    }

    console.log("TrackWrapper.render");
    return (
        <div className="cistrome-hgw-track-wrapper">
            {options.rowInfoPosition !== "hidden" ? 
                (<TrackRowInfo 
                    rowInfo={rowInfo}
                    x0={x0}
                    x1={x1}
                    y1={y1}
                    height={height}
                    infoAttrPrimary={"r1"}
                    infoAttrSecondary={"r2"}
                    rowInfoPosition={options.rowInfoPosition}
                />) : null}
        </div>
    );
}