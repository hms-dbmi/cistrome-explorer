import React, { useContext } from 'react';
import range from 'lodash/range';

import { InfoContext } from "./utils/contexts.js";
import TrackColTools from './TrackColTools.js';
import TrackRowInfo from './TrackRowInfo.js';
import TrackRowHighlight from './TrackRowHighlight.js';

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
 * @prop {(number[]|null)} selectedRows Array of row indices for selected rows. Null if all rows should be selected.
 * @prop {(number[]|null)} highlitRows Array of row indices for highlighted rows. Null if no rows should be highlighted.
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval.
 *                                          Passed down to the `TrackColTools` component.
 * @prop {function} onSelectRowInterval The function to call upon selection of a row interval.
 * @prop {function} register The function for child components to call to register their draw functions.
 */
export default function TrackWrapper(props) {
    const { 
        options, 
        multivecTrack,
        combinedTrack,
        siblingTracks,
        selectedRows,
        highlitRows,
        onSelectGenomicInterval,
        register
    } = props;

    const infoContext = useContext(InfoContext);

    if(!multivecTrack || !multivecTrack.tilesetInfo || !multivecTrack.tilesetInfo.shape) {
        // The track or track tileset info has not yet loaded.
        return null;
    }
    
    // Attributes to visualize based on the position
    const leftAttrs = options.rowInfoAttributes.filter(d => d.position === "left");
    const rightAttrs = options.rowInfoAttributes.filter(d => d.position === "right");

    const trackX = multivecTrack.position[0];
    const trackY = multivecTrack.position[1];
    const trackWidth = multivecTrack.dimensions[0];
    const trackHeight = multivecTrack.dimensions[1];
    const totalNumRows = multivecTrack.tilesetInfo.shape[1];

    // Attempt to obtain metadata values from the `tilesetInfo` field of the track.
    let rowInfo = [];
    let trackAssembly = null;
    try {
        trackAssembly = multivecTrack.tilesetInfo.coordSystem;
        // TODO: uncomment the below line to use the real metadata coming from the HiGlass Server.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        // rowInfo = multivecTrack.tilesetInfo.row_infos.map(JSON.parse);

        // TODO: remove the below line.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        rowInfo = fakedata[multivecTrack.id].tilesetInfo.rowInfo.slice(0, totalNumRows);
        if(!infoContext.state[multivecTrack.id]) {
            infoContext.dispatch({
                type: 'set_row_info',
                tilesetUid: multivecTrack.id,
                rowInfo: rowInfo
            });
        }
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
            const { field, type, order, title } = d;
            if(type === "tree") {
                // Do nothing for the "tree" type.
            } else if(type === "quantitative") {
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
            {leftAttrs.length !== 0 ? 
                (<TrackRowInfo 
                    rowInfo={transformedRowInfo}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={leftAttrs}
                    rowSort={options.rowSort}
                    rowInfoPosition="left"
                    register={register}
                />) : null}
            {rightAttrs.length !== 0 ? 
                (<TrackRowInfo
                    rowInfo={transformedRowInfo}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={rightAttrs}
                    rowSort={options.rowSort}
                    rowInfoPosition="right"
                    register={register}
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
            <TrackRowHighlight 
                trackX={trackX}
                trackY={trackY}
                trackHeight={trackHeight}
                trackWidth={trackWidth}
                totalNumRows={totalNumRows}
                selectedRows={selectedRows}
                highlitRows={highlitRows}
                register={register}
            />
        </div>
    );
}