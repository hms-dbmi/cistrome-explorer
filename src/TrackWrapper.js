import React, { useContext, useEffect, useState } from 'react';

import { InfoContext, ACTION } from "./utils/contexts.js";
import TrackRowInfo from './TrackRowInfo.js';
import TrackRowHighlight from './TrackRowHighlight.js';
import TrackRowZoomOverlay from './TrackRowZoomOverlay.js';

// TODO: remove the below fakedata import.
//       see https://github.com/hms-dbmi/cistrome-explorer/issues/26
import fakedata from './demo/fakedata/index.js';
import { getAggregatedRowInfo } from './utils/select-rows.js';

/**
 * Wrapper component associated with a particular HiGlass track.
 * @prop {object} options Options associated with the track. Contains values for all possible options.
 * @prop {object} multivecTrack A `horizontal-multivec` track object returned by `hgc.api.getTrackObject()`.
 * @prop {string} multivecTrackViewId The viewId for the multivecTrack.
 * @prop {string} multivecTrackTrackId The trackId for the multivecTrack.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a highlight interaction.
 * @prop {function} onZoomRows The function to call upon a vertical zoom interaction.
 * @prop {function} onFilterRows The function to call upon a filer interaction.
 * @prop {function} onMetadataLoad The function to call upon rowInfo is set to Context.
 * @prop {boolean} isWheelListening Whether or not to listen for wheel events for vertical zooming.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackWrapper(props) {
    const {
        options, 
        multivecTrack,
        multivecTrackViewId,
        multivecTrackTrackId,
        onAddTrack,
        onSortRows,
        onHighlightRows,
        onZoomRows,
        onFilterRows,
        onMetadataLoad,
        isWheelListening,
        drawRegister
    } = props;

    const context = useContext(InfoContext);

    const [shouldCallOnMetadataLoad, setShouldCallOnMetadataLoad] = useState(false);

    useEffect(() => {
        if(shouldCallOnMetadataLoad) {
            onMetadataLoad();
        }
    }, [shouldCallOnMetadataLoad, multivecTrackViewId, multivecTrackTrackId]);
    

    // All hooks must be above this return statement, since they need to be executed in the same order.
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
    try {
        if(["meeting-2020-04-29-track"].includes(multivecTrackTrackId)) {
            // TODO: use the below line to use the real metadata coming from the HiGlass Server tileset_info.
            //       see https://github.com/hms-dbmi/cistrome-explorer/issues/26
            rowInfo = multivecTrack.tilesetInfo.row_infos.map(JSON.parse);
        } else {
            // TODO: remove this else clause.
            //       see https://github.com/hms-dbmi/cistrome-explorer/issues/26
            rowInfo = fakedata[multivecTrack.id].tilesetInfo.rowInfo.slice(0, totalNumRows);
        }
        
        if(!context.state[multivecTrackViewId] || !context.state[multivecTrackViewId][multivecTrackTrackId]) {
            context.dispatch({
                type: ACTION.SET_ROW_INFO,
                viewId: multivecTrackViewId,
                trackId: multivecTrackTrackId,
                rowInfo: rowInfo
            });
            setShouldCallOnMetadataLoad(true);
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

    // Transformed `rowInfo` after aggregating, filtering, and sorting rows.
    // Each element of `transformedRowInfo` is either a JSON Object or an array of JSON Object
    // containing information about a single row or multiple rows that are aggregated together, respectively.
    const transformedRowInfo = (!selectedRows ? rowInfo : selectedRows.map(
        indexOrIndices => Array.isArray(indexOrIndices)
            ? rowInfo.filter((d, i) => indexOrIndices.includes(i))
            : rowInfo[indexOrIndices]
    ));
    
    // Aggregated, but not filtered, `rowInfo`.
    // This is being used for filtering interfaces since we want to allow users to filter data
    // based on the aggregated rows, not on the original rows.
    const aggregatedRowInfo = getAggregatedRowInfo(rowInfo, options.rowAggregate).map(d => d[1]);

    // console.log("TrackWrapper.render");
    return (
        <div className="cistrome-hgw-track-wrapper">
            {leftAttrs.length !== 0 ? 
                (<TrackRowInfo
                    transformedRowInfo={transformedRowInfo}
                    aggregatedRowInfo={aggregatedRowInfo}
                    viewId={multivecTrackViewId}
                    trackId={multivecTrackTrackId}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={leftAttrs}
                    rowSort={options.rowSort}
                    rowFilter={options.rowFilter}
                    rowHighlight={options.rowHighlight}
                    rowInfoPosition="left"
                    onAddTrack={onAddTrack}
                    onSortRows={onSortRows}
                    onHighlightRows={onHighlightRows}
                    onFilterRows={onFilterRows}
                    drawRegister={(key, draw, options) => {
                        drawRegister(`${key}-left`, draw, options);
                    }}
                />) : null}
            {rightAttrs.length !== 0 ? 
                (<TrackRowInfo
                    transformedRowInfo={transformedRowInfo}
                    aggregatedRowInfo={aggregatedRowInfo}
                    viewId={multivecTrackViewId}
                    trackId={multivecTrackTrackId}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={rightAttrs}
                    rowSort={options.rowSort}
                    rowFilter={options.rowFilter}
                    rowHighlight={options.rowHighlight}
                    rowInfoPosition="right"
                    onAddTrack={onAddTrack}
                    onSortRows={onSortRows}
                    onHighlightRows={onHighlightRows}
                    onFilterRows={onFilterRows}
                    drawRegister={(key, draw, options) => {
                        drawRegister(`${key}-right`, draw, options);
                    }}
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
            <TrackRowZoomOverlay
                trackX={trackX}
                trackY={trackY}
                trackHeight={trackHeight}
                trackWidth={trackWidth}
                isWheelListening={isWheelListening}
                onZoomRows={onZoomRows}
            />
        </div>
    );
}