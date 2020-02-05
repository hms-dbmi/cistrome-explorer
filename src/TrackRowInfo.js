import React, { useEffect, useRef } from 'react';
import { schemeCategory10 as d3_schemeCategory10 } from "d3-scale-chromatic";
import range from 'lodash/range';
import { mouse as d3_mouse, event as d3_event } from 'd3-selection';
import { scaleOrdinal as d3_scaleOrdinal, scaleThreshold as d3_scaleThreshold } from 'd3-scale';

import { EVENT } from './constants.js';
import { vega_scaleBand } from './utils-scales.js';
import { setupCanvas, teardownCanvas } from './utils-canvas.js';

import './TrackRowInfo.scss';

function destroyTooltip() {
    PubSub.publish(EVENT.TOOLTIP, {
        x: null,
        y: null,
        content: null
    });
}

/**
 * Component for visualization of two row info attribute values.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {string} infoAttrPrimary
 * @prop {string} infoAttrSecondary
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 */
export default function TrackRowInfo(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight, 
        rowInfo, 
        infoAttrPrimary, infoAttrSecondary,
        rowInfoPosition
    } = props;

    // Dimensions
    const top = trackY;
    const width = 120;
    const colWidth = 15;
    const xMargin = 60;
    const xMarginInitial = 5;
    const height = trackHeight;

    // Scales
    const xScale = d3_scaleThreshold();
    const yScale = vega_scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);

    const valueScalePrimary = d3_scaleOrdinal()
        .domain(Array.from(new Set(rowInfo.map(d => d[infoAttrPrimary]))))
        .range(d3_schemeCategory10);
    
    const valueScaleSecondary = d3_scaleOrdinal()
        .domain(Array.from(new Set(rowInfo.map(d => d[infoAttrSecondary]))))
        .range(d3_schemeCategory10);

    // const colorScale = d3_interpolateViridis;

    // Left offsets condition on left vs. right positioning:
    let left, primaryLeft, secondaryLeft;
    if(rowInfoPosition === "left") {
        left = trackX - xMarginInitial - width;
        primaryLeft = width - colWidth;
        secondaryLeft = width - colWidth - xMargin - colWidth;

        xScale
            .domain([secondaryLeft, secondaryLeft+colWidth, primaryLeft, primaryLeft+colWidth, width])
            .range(["m-0", "secondary", "m-1", "primary", "m-2"]);
    } else if(rowInfoPosition === "right") {
        left = trackWidth + trackX + xMarginInitial;
        primaryLeft = 0;
        secondaryLeft = colWidth + xMargin;

        xScale
            .domain([primaryLeft, primaryLeft+colWidth, secondaryLeft, secondaryLeft+colWidth, width])
            .range(["m-0", "primary", "m-1", "secondary", "m-2"]);
    }
    
    // Canvas
    const canvasRef = useRef();
    useEffect(() => {
        const { canvas, context, canvasSelection } = setupCanvas(canvasRef);
        context.clearRect(0, 0, width, height);

        // Draw primary metadata value colors.
        rowInfo.forEach((d, i) => {
            context.fillStyle = valueScalePrimary(d[infoAttrPrimary]);
            context.fillRect(primaryLeft, yScale(i), colWidth, yScale.bandwidth())
        });

        // Draw secondary metadata value colors.
        rowInfo.forEach((d, i) => {
            context.fillStyle = valueScaleSecondary(d[infoAttrSecondary]);
            context.fillRect(secondaryLeft, yScale(i), colWidth, yScale.bandwidth())
        });

        canvasSelection.on("mousemove", () => {
            const mouse = d3_mouse(canvas);
            const mouseX = mouse[0];
            const mouseY = mouse[1];

            const y = yScale.invert(mouseY);
            const x = xScale(mouseX);
            
            let xVal;
            if(x === "primary") {
                xVal = rowInfo[y][infoAttrPrimary];
            } else if(x === "secondary") {
                xVal = rowInfo[y][infoAttrSecondary];
            } else {
                destroyTooltip();
                return;
            }

            const mouseViewportX = d3_event.clientX;
            const mouseViewportY = d3_event.clientY;
            
            PubSub.publish(EVENT.TOOLTIP, {
                x: mouseViewportX,
                y: mouseViewportY,
                content: `Value: ${xVal}`
            });
        });

        canvasSelection.on("mouseout", destroyTooltip)

        return (() => teardownCanvas(canvasRef));
    });

    return (
        <div 
            className="cistrome-hgw-child"
            style={{
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height}px`,
            }}
        >
            <canvas
                ref={canvasRef}
                style={{
                    position: 'relative',
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                }}
            />
        </div>
    );
}