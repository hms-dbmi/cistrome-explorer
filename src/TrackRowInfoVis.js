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
 * @prop {object} rowSort The options for sorting rows.
 * @prop {object} rowFilter The options for filtering rows.
 * @prop {object} rowHighlight The options for highlighting rows.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 * @prop {function} onWidthChanged The function to call when the component resize element has been dragged.
 */
export default function TrackRowInfoVis(props) {
    const {
        left, top, width, height,
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


    const divRef = useRef();
    const resizerRef = useRef();

    const [isHovering, setIsHovering] = useState(false);
    const dragX = useRef(null);

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const event = d3.event;
        dragX.current = event.sourceEvent.clientX;
    }, [dragX])

    const ended = useCallback(() => {
        dragX.current = null;
    }, [dragX])

    const dragged = useCallback(() => {
        const event = d3.event;
        const diff = event.sourceEvent.clientX - dragX.current;
        let newWidth = isLeft ? width - diff : width + diff;
        if(newWidth < minWidth) {
            newWidth = minWidth;
        }
        // Emit the new width value to the parent component.
        onWidthChanged(newWidth);
    }, [dragX, width, onWidthChanged]);

    // Detect drag events for the resize element.
    useEffect(() => {
        const resizer = resizerRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(resizer).call(drag);

        return () => d3.select(resizer).on(".drag", null);
    }, [resizerRef, started, dragged, ended]);

    // Detect hover events to know when to show the resize element.
    useEffect(() => {
        const div = divRef.current;

        d3.select(div).on("mouseenter", () => setIsHovering(true));
        d3.select(div).on("mouseleave", () => setIsHovering(false));

        return () => {
            d3.select(div).on("mouseenter", null);
            d3.select(div).on("mouseleave", null);
        };
    }, [divRef, setIsHovering]);

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

    // Create the resizer element.
    const resizer = useMemo(() => {
        const lrKey = (isLeft ? "left" : "right");
        return (
            <div
                ref={resizerRef}
                className="visualization-resizer"
                style={{
                    top: `${top + (height - resizerHeight) / 2.0}px`,
                    [lrKey]: resizerMargin,
                    height: `${resizerHeight}px`,
                    width: `${resizerWidth}px`,
                    opacity: (isHovering ? 1 : 0)
                }}
            />
        );
    }, [isLeft, top, height, isHovering]);

    return (
        <div 
            ref={divRef}
            style={{
                position: 'absolute',
                top: `${top}px`,
                left: `${left}px`, 
                height: `${height}px`,
            }}
        >
            {React.createElement(
                fieldTypeToVisComponent[fieldInfo.type],
                {
                    left,
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