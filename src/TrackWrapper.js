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
 * @prop {(object|null)} combinedTrack A `combined` track object returned by `hgc.api.getTrackObject()`. 
 *                              If not null, it is the parent track of the `multivecTrack`.
 * @prop {object[]} siblingTracks An array of `viewport-projection-horizontal` track objects, which
 *                                are siblings of `multivecTrack` (children of the same `combined` track).
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval. 
 *                                          Passed down to the `TrackColTools` component.
 * @prop {function} register The function for child components to call to register their draw functions.
 */
export default function TrackWrapper(props) {
    const { 
        options, 
        multivecTrack,
        combinedTrack,
        siblingTracks,
        onSelectGenomicInterval,
        register
    } = props;

    if(!multivecTrack || !multivecTrack.tilesetInfo) {
        // The track or track tileset info has not yet loaded.
        return null;
    }

    const trackX = multivecTrack.position[0];
    const trackY = multivecTrack.position[1];
    const trackWidth = multivecTrack.dimensions[0];
    const trackHeight = multivecTrack.dimensions[1];

    // Attempt to obtain metadata values from the `tilesetInfo` field of the track.
    let rowInfo = [];
    let trackAssembly = null;
    try {
        trackAssembly = multivecTrack.tilesetInfo.coordSystem;
        // TODO: uncomment the below line to use the real metadata coming from the HiGlass Server.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        // rowInfo = multivecTrack.tilesetInfo.row_infos.map(JSON.parse);

        // TODO: remove the below lines.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        const numRows = multivecTrack.tilesetInfo.shape[1];
        rowInfo = fakedata[multivecTrack.id].tilesetInfo.rowInfo.slice(0, numRows);
    } catch(e) {
        console.log(e);
    }

    /*
     * Transform data based on options (e.g., sorting, filtering).
     */
    // Filter
    // ...
    
    // Sort
    let transformedRowInfo = rowInfo.slice();
    if(options.rowSort && options.rowSort.length > 0) {
        let sortOptions = options.rowSort.slice().reverse();
        sortOptions.forEach((d, i) => {
            const { field, type, order } = d;
            if(type === "quantitative") {
                transformedRowInfo.sort((a, b) => (a[field] - b[field]) * (order === "ascending" ? 1 : -1));
            } else {
                transformedRowInfo.sort(function(a, b) {
                    let compared = 0, categoryA = a[field].toUpperCase(), categoryB = b[field].toUpperCase();
                    if(categoryA > categoryB) {
                        compared = 1;
                    } else {
                        compared = -1;
                    }
                    return compared * (order === "ascending" ? 1 : -1);
                });
            }
        });
    }

    console.log("TrackWrapper.render");
    return (
        <div className="cistrome-hgw-track-wrapper">
            {options.rowInfoPosition !== "hidden" ? 
                (<TrackRowInfo 
                    rowInfo={transformedRowInfo}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={options.rowInfoAttributes}
                    rowInfoPosition={options.rowInfoPosition}
                    rowSort={options.rowSort}
                    register={register}
                />) : null}
            {options.rowLinkPosition !== "hidden" ? 
                (<TrackRowLink
                    rowInfo={transformedRowInfo}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowLinkAttribute={options.rowLinkAttribute}
                    rowLinkNameAttribute={options.rowLinkNameAttribute}
                    rowLinkPosition={options.rowLinkPosition}
                />) : null}
            {options.colToolsPosition !== "hidden" ? 
                (<TrackColTools
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    trackAssembly={trackAssembly}
                    combinedTrack={combinedTrack}
                    siblingTracks={siblingTracks}
                    colToolsPosition={options.colToolsPosition}
                    onSelectGenomicInterval={onSelectGenomicInterval}
                    register={register}
                />) : null}
        </div>
    );
}