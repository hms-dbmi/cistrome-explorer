import React, { useState, useCallback } from 'react';
import d3 from './utils/d3.js';

import TrackRowInfoVis from "./TrackRowInfoVis.js";
import { selectRows } from './utils/select-rows.js';

const DEFAULT_UNIT_WIDTH = 100;
const DEFAULT_BAND_WIDTH = 60;

/**
 * Parent component for visualization of multiple row info attribute values, on a particular side (left/right).
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
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
        originalRowInfo,    // TODO: Better way to differentiate the tree
        rowInfo,
        transformedRowInfo, 
        selectedRows, // TODO: 
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

    // Dimensions
    const isLeft = rowInfoPosition === "left";

    const top = trackY;
    const height = trackHeight;
    const [unitWidths, setUnitWidths] = useState(rowInfoAttributes.map(() => DEFAULT_UNIT_WIDTH));
    
    const setTrackProps = useCallback(() => {
        const properties = [];
        
        let currentLeft = 0, prevTrackSelectedRows = selectedRows;
        let isPrevIndependentYScale = false;
        rowInfoAttributes.forEach((attribute, i) => {
            const fieldInfo = isLeft ? rowInfoAttributes[rowInfoAttributes.length - i - 1] : attribute;
            const { field, type, resolveYScale, sort: order } = fieldInfo;        
            const width = unitWidths[i];
            let isCurrIndependentYScale = false;
            
            // Determine whether to use a track-level transformedRowInfo or the global one
            let currSelectedRows, currTransformedRowInfo;
            if(resolveYScale && order) {
                isCurrIndependentYScale = true;

                // Track-independent arrangement
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
            
            if(isPrevIndependentYScale || isCurrIndependentYScale) {
                // We need to show band connections before this track.
                properties.push({
                    top: 0,
                    left: currentLeft, 
                    width: DEFAULT_BAND_WIDTH,
                    height,
                    fieldInfo: { type: "band" },
                    leftSelectedRows: prevTrackSelectedRows,
                    rightSelectedRows: currSelectedRows
                });
                currentLeft += DEFAULT_BAND_WIDTH;
            }

            // Add props for the current vertical track
            properties.push({
                index: i,
                top: 0,
                left: currentLeft, 
                width,
                height,
                fieldInfo,
                transformedRowInfo: currTransformedRowInfo
            });
            currentLeft += width;

            prevTrackSelectedRows = currSelectedRows;
            isPrevIndependentYScale = isCurrIndependentYScale;
        });

        if(isLeft && isPrevIndependentYScale) {
            // Add a last band track between a HiGlass heatmap and a vertical track on the left
            properties.push({
                top: 0,
                left: currentLeft, 
                width: DEFAULT_BAND_WIDTH,
                height,
                fieldInfo: { type: "band" },
                leftSelectedRows: prevTrackSelectedRows,
                rightSelectedRows: selectedRows
            });
        }
        return properties;
    }, [unitWidths, rowInfoAttributes, selectedRows, originalRowInfo, transformedRowInfo]);

    const trackProps = setTrackProps();
    const totalWidth = d3.sum(trackProps.map(d => d.width));
    const left = isLeft ? trackX - totalWidth : trackX + trackWidth;

    const setUnitWidthByIndex = useCallback((i, val) => {
        console.log(i, val); // TODO:
        unitWidths[i] = val;
        setUnitWidths([...unitWidths]);
    }, [unitWidths]);

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
                    onWidthChanged={(val) => setUnitWidthByIndex(d.index, val)}
                />
            ))}
        </div>
    );
}