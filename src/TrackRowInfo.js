import React, { useState } from 'react';
import d3 from './utils/d3.js';

import TrackRowInfoVis from "./TrackRowInfoVis.js";

/**
 * Parent component for visualization of multiple row info attribute values, on a particular side (left/right).
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {object[]} aggregatedRowInfo The `rowInfo` array after aggregating rows.
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
        transformedRowInfo, 
        aggregatedRowInfo,
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
    const defaultUnitWidth = 100;
    const [unitWidths, setUnitWidths] = useState(rowInfoAttributes.map(() => defaultUnitWidth));
    
    const totalWidth = d3.sum(unitWidths);
    const left = isLeft ? trackX - totalWidth : trackX + trackWidth;

    function setUnitWidthByIndex(i, val) {
        const newUnitWidths = Array.from(unitWidths);
        newUnitWidths[i] = val;
        setUnitWidths(newUnitWidths);
    }
    
    // Determine position of each dimension.
    let trackProps = [], currentLeft = 0;
    rowInfoAttributes.forEach((attribute, i) => {
        const fieldInfo = isLeft ? rowInfoAttributes[rowInfoAttributes.length - i - 1] : attribute;
        const width = unitWidths[i];
        trackProps.push({
            top: 0,
            left: currentLeft, 
            width,
            height,
            fieldInfo,
        });
        currentLeft += width;
    });

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
                    transformedRowInfo={transformedRowInfo}
                    aggregatedRowInfo={aggregatedRowInfo}
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