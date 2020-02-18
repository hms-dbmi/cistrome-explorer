import React, { useContext } from 'react';
import range from 'lodash/range';

import { InfoContext, ACTION } from "./utils/contexts.js";
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
 * @prop {string} multivecTrackViewId The viewId for the multivecTrack.
 * @prop {string} multivecTrackTrackId The trackId for the multivecTrack.
 * @prop {(object|null)} combinedTrack A `combined` track object returned by `hgc.api.getTrackObject()`.
 *                              If not null, it is the parent track of the `multivecTrack`.
 * @prop {object[]} siblingTracks An array of `viewport-projection-horizontal` track objects, which
 *                                are siblings of `multivecTrack` (children of the same `combined` track).
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval.
 *                                          Passed down to the `TrackColTools` component.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackWrapper(props) {
    const { 
        options, 
        multivecTrack,
        multivecTrackViewId,
        multivecTrackTrackId,
        combinedTrack,
        siblingTracks,
        onSelectGenomicInterval,
        onSortRows,
        onSearchRows,
        drawRegister
    } = props;

    const context = useContext(InfoContext);

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
        
        if(!context.state[multivecTrackViewId] || !context.state[multivecTrackViewId][multivecTrackTrackId]) {
            context.dispatch({
                type: ACTION.SET_ROW_INFO,
                viewId: multivecTrackViewId,
                trackId: multivecTrackTrackId,
                rowInfo: rowInfo
            });
        }
    } catch(e) {
        console.log(e);
    }

    let selectedRows;
    let highlitRows;
    try {
        selectedRows = context.state[multivecTrackViewId][multivecTrackTrackId].selectedRows;
        highlitRows = context.state[multivecTrackViewId][multivecTrackTrackId].highlitRows;
    } catch(e) {
        // pass
        console.log(e);
    }

    const transformedRowInfo = (!selectedRows ? rowInfo : selectedRows.map(i => rowInfo[i]));

    console.log("TrackWrapper.render");
    return (
        <div className="cistrome-hgw-track-wrapper">
            {leftAttrs.length !== 0 ? 
                (<TrackRowInfo 
                    rowInfo={transformedRowInfo}
                    viewId={multivecTrackViewId}
                    trackId={multivecTrackTrackId}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={leftAttrs}
                    rowSort={options.rowSort}
                    rowInfoPosition="left"
                    onSortRows={onSortRows}
                    onSearchRows={onSearchRows}
                    drawRegister={drawRegister}
                />) : null}
            {rightAttrs.length !== 0 ? 
                (<TrackRowInfo
                    rowInfo={transformedRowInfo}
                    viewId={multivecTrackViewId}
                    trackId={multivecTrackTrackId}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={rightAttrs}
                    rowSort={options.rowSort}
                    rowInfoPosition="right"
                    onSortRows={onSortRows}
                    onSearchRows={onSearchRows}
                    drawRegister={drawRegister}
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
                    drawRegister={drawRegister}
                />) : null}
            <TrackRowHighlight 
                trackX={trackX}
                trackY={trackY}
                trackHeight={trackHeight}
                trackWidth={trackWidth}
                totalNumRows={totalNumRows}
                selectedRows={selectedRows}
                highlitRows={highlitRows}
                drawRegister={drawRegister}
            />
        </div>
    );
}