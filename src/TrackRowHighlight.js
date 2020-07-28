import React, { useEffect, useRef, useCallback } from 'react';
import range from 'lodash/range';
import d3 from './utils/d3.js';
import Two from './utils/two.js';
import { HIGHLIGHTING_STROKE, HIGHLIGHTING_COLOR, HIGHLIGHTING_OPACITY } from './utils/linking-views.js';

/**
 * Component for visualizing highlighted rows of a multivec track.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {number} totalNumRows The number of rows.
 * @prop {(number[]|null)} selectedRows Array of selected row indices.
 * @prop {(number[]|null)} highlitRows Array of highlighted row indices.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowHighlight(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight,
        totalNumRows,
        selectedRows,
        highlitRows,
        drawRegister
    } = props;

    const top = trackY;
    const left = trackX;
    const height = trackHeight;
    const width = trackWidth;

    // Render canvas
    const canvasRef = useRef();

    const yScale = d3.scaleBand()
        .range([0, height]);

    if(selectedRows) {
        // Non-null selectedRows means a subset of rows are selected, potentially in a particular order.
        yScale.domain(selectedRows);
    } else {
        // Null selectedRows means all rows are selected, in the default order.
        yScale.domain(range(totalNumRows));
    }

    const rowHeight = yScale.bandwidth();

    // Render each track.
    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        if(highlitRows) {
            let aggregatedRows = 1;
            const sortedHighlitRows = highlitRows.slice().sort((a, b) => selectedRows.indexOf(a) - selectedRows.indexOf(b));
            sortedHighlitRows.forEach((d, i) => {
                const startY = yScale(d);
                
                if(
                    sortedHighlitRows.length > i + 1 && 
                    Math.abs(selectedRows.indexOf(sortedHighlitRows[i + 1]) - selectedRows.indexOf(d)) === 1
                ) {
                    // Aggregate highlighting rows for drawing the stroke only once
                    aggregatedRows++;
                    return;
                }

                const rect = two.makeRect(0, startY - rowHeight * (aggregatedRows - 1), width, rowHeight * aggregatedRows);
                rect.stroke = HIGHLIGHTING_STROKE;
                rect.fill = HIGHLIGHTING_COLOR;
                rect.opacity = HIGHLIGHTING_OPACITY;

                aggregatedRows = 1;
            });
        }

        two.update();
        return two.teardown;
    }, [width, height, totalNumRows, selectedRows, highlitRows]);

    drawRegister("TrackRowHighlight", draw, { top, left, width, height });

    useEffect(() => {
        const canvas = canvasRef.current;
        const teardown = draw(canvas);
        return teardown;
    });

    //console.log("TrackRowHighlight.render");
    return (
        <canvas
            ref={canvasRef}
            style={{
                position: 'absolute',
                top: `${top}px`,
                left: `${left}px`,
                width: `${width}px`,
                height: `${height}px`,
                pointerEvents: 'none'
            }}
        />
    );
}