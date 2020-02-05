import React from 'react';

import TrackColTools from './TrackColTools.js';
import TrackRowInfo from './TrackRowInfo.js';
import TrackRowLink from './TrackRowLink.js';

// TODO: remove the below fakedata import.
//       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
import fakedata from './demo/fakedata/index.js';

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
        // TODO: uncomment the below line to use the real metadata coming from the HiGlass Server.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        // rowInfo = multivecTrack.tilesetInfo.row_infos.map(JSON.parse);

        // TODO: remove the below lines.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        const numRows = multivecTrack.tilesetInfo.shape[1];
        rowInfo = fakedata[multivecTrack.id].tilesetInfo.rowInfo.slice(0, numRows);
    } catch(e) {
        console.log(e)
    }

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
                    infoAttributes={options.infoAttributes}
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
                    rowLinkAttribute={options.rowLinkNameAttribute}
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