import React from 'react';

import TrackColTools from './TrackColTools.js';
import TrackRowInfo from './TrackRowInfo.js';
import TrackRowLink from './TrackRowLink.js';

/**
 * Wrapper component associated with a particular HiGlass track.
 * @prop {object} options Options associated with the track. Contains values for all possible options.
 * @prop {object} track A track object returned by `hgc.api.getTrackObject()`.
 */
export default function TrackWrapper(props) {
    const { 
        options, 
        track,
        onSelectGenomicInterval
    } = props;

    if(!track || !track.tilesetInfo) {
        // The track or track tileset info has not yet loaded.
        return null;
    }

    const trackX = track.position[0];
    const trackY = track.position[1];
    const trackWidth = track.dimensions[0];
    const trackHeight = track.dimensions[1];
    let rowInfo = [];
    try {
        rowInfo = track.tilesetInfo.row_infos.map(JSON.parse);
    } catch(e) {
        // TODO: Remove this catch block. This is only being used for the resgen.io cistrome demo data, 
        //       since its metadata is stored in tab-separated strings rather than JSON objects.
        rowInfo = track.tilesetInfo.row_infos.map(d => d.split("\t"));
    }

    // TODO: Obtain infoAttr keys from `options`.

    console.log("TrackWrapper.render");
    return (
        <div className="cistrome-hgw-track-wrapper">
            {options.rowInfoPosition !== "hidden" ? 
                (<TrackRowInfo 
                    rowInfo={rowInfo}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    infoAttrPrimary={4}
                    infoAttrSecondary={1}
                    rowInfoPosition={options.rowInfoPosition}
                />) : null}
            {options.rowLinkPosition !== "hidden" ? 
                (<TrackRowLink
                    rowInfo={rowInfo}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowLinkAttribute={options.rowLinkAttribute}
                    rowLinkPosition={options.rowLinkPosition}
                />) : null}
            {options.colToolsPosition !== "hidden" ? 
                (<TrackColTools
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    colToolsPosition={options.colToolsPosition}
                    onSelectGenomicInterval={onSelectGenomicInterval}
                />) : null}
        </div>
    );
}