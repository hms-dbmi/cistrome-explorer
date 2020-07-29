import React, { useState, useEffect, useCallback } from 'react';
import d3 from './utils/d3.js';

import TrackRowInfoVis from "./TrackRowInfoVis.js";
import { selectRows } from './utils/select-rows.js';
import { modifyItemInArray } from './utils/array.js';
import { getNumOfTracks } from './utils/layout.js';

const DEFAULT_TRACK_WIDTH = 100;
const DEFAULT_BAND_WIDTH = 100;

const TRACK_TYPE = Object.freeze({
    VIS: "vis",
    BAND: "band"
});

/**
 * Parent component for visualization of multiple row info attribute values, on a particular side (left/right).
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {object[]} originalRowInfo The array of JSON Object containing `raw` row information.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {array} selectedRows The array of selected indices. 
 * @prop {array} highlitRows The array of highlit indices.
 * @prop {array} rowInfoAttributes Array of JSON object, one object for the names and types of each attribute.
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 * @prop {object} rowSort The options for sorting rows.
 * @prop {object} rowFilter The options for filtering rows.
 * @prop {object} rowHighlight The options for highlighting rows.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a highlight interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfo(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight,
        originalRowInfo, // Being used only for determining track-specific `selectedRows`
        rowInfo, // Being used only for showing data in filtering interfaces (e.g., keyword search results)
        transformedRowInfo, // Being used for visualizations
        selectedRows,
        highlitRows,
        rowInfoAttributes,
        rowInfoPosition,
        rowSort,
        rowFilter,
        rowHighlight,
        onAddTrack,
        onSortRows,
        onHighlightRows,
        onFilterRows,
        drawRegister
    } = props;

    const isLeft = rowInfoPosition === "left";
    const top = trackY;
    const height = trackHeight;
    const [trackWidths, setTrackWidths] = useState(new Array(
        getNumOfTracks(isLeft, rowInfoAttributes.map(d => d.resolveYScale))
    ));
    const trackProps = generateTrackProps();

    // This function generates properties for tracks (e.g., width, starting position),
    // considering the band-connection tracks between actual tracks
    function generateTrackProps() { 
        const properties = [];
        
        let currentLeft = 0, prevTrackSelectedRows = selectedRows;
        let isPrevIndependentYScale = false;
        let currentIndex = 0;

        rowInfoAttributes.forEach((attribute, i) => {
            // select the `fieldInfo` that is closest to HiGlass first
            const fieldInfo = isLeft ? rowInfoAttributes[rowInfoAttributes.length - i - 1] : attribute;
            const { field, type, resolveYScale, sort: order, width: initWidth } = fieldInfo;
            let isCurrIndependentYScale = false;
            
            // Determine whether to use a track-specific transformedRowInfo or the global one
            let currSelectedRows, currTransformedRowInfo;
            if(resolveYScale && order) {
                isCurrIndependentYScale = true;

                // Use track-specific arrangement
                currSelectedRows = selectRows(
                    originalRowInfo,
                    { rowFilter },
                    { rowSort: [{ field, type, order }] }
                );
                currTransformedRowInfo = (!currSelectedRows ? originalRowInfo : currSelectedRows.map(
                    indexOrIndices => Array.isArray(indexOrIndices)
                        ? originalRowInfo.filter((d, i) => indexOrIndices.includes(i))
                        : originalRowInfo[indexOrIndices]
                ));
            } else {
                currSelectedRows = selectedRows;
                currTransformedRowInfo = transformedRowInfo;
            }
            
            if(
                (isPrevIndependentYScale || isCurrIndependentYScale)
                && !(isLeft && i === 0) // we don't want to add a band-connection track on the left-most area
            ) {
                // We need to show band connections before this regular track
                const width = trackWidths[currentIndex++] ?? DEFAULT_BAND_WIDTH;
                properties.push({
                    type: TRACK_TYPE.BAND,
                    top: 0,
                    left: currentLeft, 
                    width,
                    height,
                    fieldInfo: { type: "band" },
                    leftSelectedRows: prevTrackSelectedRows,
                    rightSelectedRows: currSelectedRows
                });
                currentLeft += width;
            }

            // Add props for the current regular track
            const width = trackWidths[currentIndex++] ?? initWidth ?? DEFAULT_TRACK_WIDTH;
            properties.push({
                type: TRACK_TYPE.VIS,
                top: 0,
                left: currentLeft, 
                width,
                height,
                fieldInfo,
                selectedRows: currSelectedRows,
                transformedRowInfo: currTransformedRowInfo
            });
            currentLeft += width;

            prevTrackSelectedRows = currSelectedRows;
            isPrevIndependentYScale = isCurrIndependentYScale;
        });

        if(isLeft && isPrevIndependentYScale) {
            // Add a last band track between a HiGlass heatmap and a vis track on the left
            properties.push({
                type: TRACK_TYPE.BAND,
                top: 0,
                left: currentLeft, 
                width: trackWidths[currentIndex] ?? DEFAULT_BAND_WIDTH,
                height,
                fieldInfo: { type: "band" },
                leftSelectedRows: prevTrackSelectedRows,
                rightSelectedRows: selectedRows
            });
        }
        return properties;
    }

    const totalWidth = d3.sum(trackProps.map(d => d.width));
    const left = isLeft ? trackX - totalWidth : trackX + trackWidth;

    const setUnitWidthByIndex = useCallback((i, val) => {
        const newTrackWidths = modifyItemInArray(trackWidths, i, val);
        setTrackWidths(newTrackWidths);
    }, [trackWidths]);

    //console.log("TrackRowInfo.render");
    return (
        <div
            style={{
                position: 'absolute',
                top: `${top}px`,
                left: `${left}px`,
                height: `${height}px`
            }}
        >
            {trackProps.map((d, i) => (
                <TrackRowInfoVis 
                    key={i}
                    top={d.top}
                    left={d.left}
                    width={d.width}
                    height={d.height}
                    isLeft={isLeft}
                    fieldInfo={d.fieldInfo}
                    rowInfo={rowInfo}
                    transformedRowInfo={d.transformedRowInfo}
                    leftSelectedRows={d.leftSelectedRows}
                    rightSelectedRows={d.rightSelectedRows}
                    selectedRows={d.selectedRows}
                    highlitRows={highlitRows}
                    rowSort={rowSort}
                    rowFilter={rowFilter}
                    rowHighlight={rowHighlight}
                    onAddTrack={onAddTrack}
                    onSortRows={onSortRows}
                    onHighlightRows={onHighlightRows}
                    onFilterRows={onFilterRows}
                    drawRegister={(key, draw, options) => {
                        drawRegister(`${key}-${i}`, draw, { top: top + d.top, left: left + d.left, width: d.width, height: d.height })
                    }}
                    onWidthChanged={(val) => setUnitWidthByIndex(i, val)}
                />
            ))}
        </div>
    );
}