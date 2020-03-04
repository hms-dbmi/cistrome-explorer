import React, { useState, useEffect, useCallback, useMemo } from 'react';
import d3 from './utils/d3.js';
import { modifyItemInArray } from './utils/array.js';

import TrackRowInfoVis from "./TrackRowInfoVis.js";

/**
 * Parent component for visualization of multiple row info attribute values, on a particular side (left/right).
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {array} rowInfoAttributes Array of JSON object, one object for the names and types of each attribute.
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 * @prop {object} rowSort The options for sorting rows.
 * @prop {object} rowFilter The options for filtering rows.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfo(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight, 
        transformedRowInfo, 
        rowInfoAttributes,
        rowInfoPosition,
        rowSort,
        rowFilter,
        rowHighlight,
        onSortRows,
        onSearchRows,
        onFilterRows,
        drawRegister
    } = props;

    // Dimensions
    const isLeft = rowInfoPosition === "left";

    const top = trackY;
    const height = trackHeight;
    const defaultUnitWidth = 100;
    const [unitWidths, setUnitWidths] = useState(rowInfoAttributes.map(() => defaultUnitWidth));


    const lrKey = (isLeft ? "right" : "left");
    const lrVal = (isLeft ? trackX + trackWidth : trackX + trackWidth);

    function setUnitWidthByIndex(i, val) {
        const newUnitWidths = Array.from(unitWidths);
        newUnitWidths[i] = val;
        setUnitWidths(newUnitWidths);
    }

    
    
    // Determine position of each dimension.
    let trackProps = [], currentLr = 0;
    rowInfoAttributes.forEach((attribute, i) => {
        const fieldInfo = attribute;
        const width = unitWidths[i];
        trackProps.push({
            top: 0,
            [lrKey]: currentLr, 
            width,
            height,
            fieldInfo,
        });
        currentLr += width;
    });

    console.log("TrackRowInfo.render");
    return (
        <div
            style={{
                border: '1px solid green',
                position: 'absolute',
                top: `${top}px`,
                [lrKey]: `${lrVal}px`,
                height: `${height}px`
            }}
        >
            {trackProps.map((d, i) => (
                <TrackRowInfoVis 
                    key={i}
                    top={d.top}
                    left={d.left}
                    right={d.right}
                    width={d.width}
                    height={d.height}
                    isLeft={isLeft}
                    fieldInfo={d.fieldInfo}
                    transformedRowInfo={transformedRowInfo}
                    rowSort={rowSort}
                    rowFilter={rowFilter}
                    rowHighlight={rowHighlight}
                    onSortRows={onSortRows}
                    onSearchRows={onSearchRows}
                    onFilterRows={onFilterRows}
                    drawRegister={drawRegister}
                    
                    onWidthChanged={(val) => setUnitWidthByIndex(i, val)}
                />
            ))}
        </div>
    );
}