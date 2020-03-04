import React, { useRef, useState, useEffect, useCallback, useMemo } from 'react';
import d3 from './utils/d3.js';

import TrackRowInfoVisNominalBar from './TrackRowInfoVisNominalBar.js';
import TrackRowInfoVisQuantitativeBar from './TrackRowInfoVisQuantitativeBar.js';
import TrackRowInfoVisLink from './TrackRowInfoVisLink.js';
import TrackRowInfoVisDendrogram from './TrackRowInfoVisDendrogram.js';
import './TrackResizer.scss';


const fieldTypeToVisComponent = {
    "nominal": TrackRowInfoVisNominalBar,
    "quantitative": TrackRowInfoVisQuantitativeBar,
    "url": TrackRowInfoVisLink,
    "tree": TrackRowInfoVisDendrogram
};

/**
 * General component for visualization of a particular row info attribute.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} transformedRowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {string} titleSuffix The suffix of a title, information about sorting and filtering status.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVis(props) {
    const {
        left, right, top, width, height,
        fieldInfo,
        isLeft,
        transformedRowInfo,
        rowSort,
        rowFilter,
        rowHighlight,
        onSortRows,
        onSearchRows,
        onFilterRows,
        drawRegister,
        onWidthChanged
    } = props;

    const minWidth = 40;
    const resizerWidth = 4
    const resizerHeight = 10
    const resizerMargin = 2;
    const marginForMouseEvent = 20;

    const lrKey = (isLeft ? "right" : "left");
    const lrVal = (isLeft ? right : left);

    const resizerRef = useRef();

    const [isResizing, setIsResizing] = useState(false);

    const dragged = useCallback(() => {
        const event = d3.event;
        const mouseX = (isLeft ? event.x : event.x);

        let newWidth = isLeft ? width - mouseX : mouseX;
        if(newWidth < minWidth) {
            newWidth = minWidth;
        }
        console.log(event.x, event.dx)
        console.log(newWidth);
        // Emit the new width value to the parent component.
        onWidthChanged(newWidth);
    }, [onWidthChanged]);

    useEffect(() => {
        const resizer = resizerRef.current;

        console.log("effect")

        const drag = d3.drag()
            .on("start", () => setIsResizing(true))
            .on("drag", dragged)
            .on("end", () => setIsResizing(false));

        d3.select(resizer).call(drag);

    }, [dragged]);

    // Determine the title suffix.
    let titleSuffix = "";
    const sortInfo = rowSort ? rowSort.find(d => d.field === fieldInfo.field) : undefined;
    if(sortInfo) {
        titleSuffix += ` | sorted (${sortInfo.order})`;
    }
    const filterInfo = rowFilter ? rowFilter.filter(d => d.field === fieldInfo.field) : undefined;
    if(filterInfo && filterInfo.length > 0) {
        titleSuffix += ` | filtered by ${filterInfo.map(d => `"${d.contains}"`).join(', ')}`;
    }
    if(rowHighlight && rowHighlight.field === fieldInfo.field && rowHighlight.contains.length > 0) {
        titleSuffix += ` | highlighted by "${rowHighlight.contains}"`;
    }

    const resizer = (
        <div
            ref={resizerRef}
            className="visualization-resizer"
            style={{
                top: `${top + (height - resizerHeight) / 2.0}px`,
                left: `${isLeft ? resizerMargin : width - resizerWidth - resizerMargin}px`,
                height: `${resizerHeight}px`,
                width: `${resizerWidth}px`
            }}
        />
    );

    return (
        <div 
            style={{
                border: '1px solid red',
                position: 'absolute',
                top: `${top}px`,
                [lrKey]: `${lrVal}px`, 
                height: `${height}px`,
            }}
        >
            {React.createElement(
                fieldTypeToVisComponent[fieldInfo.type],
                {
                    [lrKey]: lrVal,
                    top: 0,
                    width,
                    height,
                    isLeft,
                    fieldInfo,
                    transformedRowInfo,
                    titleSuffix,
                    onSortRows,
                    onSearchRows,
                    onFilterRows,
                    drawRegister,
                }
            )}
            {resizer}
        </div>
    );
}