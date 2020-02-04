import React from 'react';

import TrackColTools from './TrackColTools.js';
import TrackRowInfo from './TrackRowInfo.js';
import TrackRowLink from './TrackRowLink.js';

/**
 * Wrapper component associated with a particular HiGlass track.
 * @prop {object} options Options associated with the track. Contains values for all possible options.
 * @prop {object} multivecTrack A `horizontal-multivec` track object returned by `hgc.api.getTrackObject()`.
 * @prop {object} combinedTrack A `combined` track object returned by `hgc.api.getTrackObject()`. 
 *                              If not null, it is the parent track of the `multivecTrack`.
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval. 
 *                                          Passed down to the `TrackColTools` component.
 */
export default function TrackWrapper(props) {
    const { 
        options, 
        multivecTrack,
        combinedTrack,
        onSelectGenomicInterval
    } = props;

    if(!multivecTrack || !multivecTrack.tilesetInfo) {
        // The track or track tileset info has not yet loaded.
        return null;
    }

    const trackX = multivecTrack.position[0];
    const trackY = multivecTrack.position[1];
    const trackWidth = multivecTrack.dimensions[0];
    const trackHeight = multivecTrack.dimensions[1];
    let rowInfo = [];
    try {
        rowInfo = multivecTrack.tilesetInfo.row_infos.map(JSON.parse);
    } catch(e) {
        // TODO: Remove this catch block. This is only being used for the resgen.io cistrome demo data, 
        //       since its metadata is stored in tab-separated strings rather than JSON objects.
        rowInfo = multivecTrack.tilesetInfo.row_infos.map(d => d.split("\t"));
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
                    infoAttrPrimary={options.infoAttrPrimary}
                    infoAttrSecondary={options.infoAttrSecondary}
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
            {(!combinedTrack && options.colToolsPosition !== "hidden") ? 
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