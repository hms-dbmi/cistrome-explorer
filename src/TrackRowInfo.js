import React, { useEffect, useRef } from 'react';
import { scale as vega_scale } from 'vega-scale';
import { interpolateViridis as d3_interpolateViridis } from "d3-scale-chromatic";
import range from 'lodash/range';
import { mouse as d3_mouse, event as d3_event } from 'd3';

import { setupCanvas, teardownCanvas } from './utils-canvas.js';

import './TrackRowInfo.scss';

/**
 * Component for visualization of two row info attribute values.
 * @prop {number} x0
 * @prop {number} x1
 * @prop {number} y0
 * @prop {number} y1
 * @prop {number} height The track height.
 * @prop {Array} rowInfo Array of JSON objects, one object for each row.
 * @prop {string} infoAttrPrimary
 * @prop {string} infoAttrSecondary
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 */
export default function TrackRowInfo(props) {

    const {
        x0, x1, y1, height, rowInfo, 
        infoAttrPrimary, infoAttrSecondary,
        rowInfoPosition
    } = props;

    const top = y1;
    const width = 120;
    const colWidth = 15;
    const xMargin = 60;
    const xMarginInitial = 5;

    // Left offsets condition on left vs. right positioning:
    let left, primaryLeft, secondaryLeft;
    if(rowInfoPosition === "left") {
        left = x1 - xMarginInitial - width;
        primaryLeft = width - colWidth;
        secondaryLeft = width - colWidth - xMargin - colWidth;
    } else if(rowInfoPosition === "right") {
        left = x0 + x1 + xMarginInitial;
        primaryLeft = 0;
        secondaryLeft = colWidth + xMargin;
    }
    

    // Scales
    const yScale = vega_scale('band')()
        .domain(range(rowInfo.length - 1))
        .range([0, height]);

    const valueScalePrimary = vega_scale('band')()
        .domain(Array.from(new Set(rowInfo.map(d => d[infoAttrPrimary]))))
        .range([0, 1]);
    
    const valueScaleSecondary = vega_scale('band')()
        .domain(Array.from(new Set(rowInfo.map(d => d[infoAttrSecondary]))))
        .range([0, 1]);

    const colorScale = d3_interpolateViridis;
    
    // Canvas
    const canvasRef = useRef();
    useEffect(() => {
        const { canvas, context, canvasSelection } = setupCanvas(canvasRef);
        context.clearRect(0, 0, width, height);

        // Draw primary metadata value colors.
        rowInfo.forEach((d, i) => {
            context.fillStyle = colorScale(valueScalePrimary(d[infoAttrPrimary]));
            context.fillRect(primaryLeft, yScale(i), colWidth, yScale.bandwidth())
        });

        // Draw secondary metadata value colors.
        rowInfo.forEach((d, i) => {
            context.fillStyle = colorScale(valueScaleSecondary(d[infoAttrSecondary]));
            context.fillRect(secondaryLeft, yScale(i), colWidth, yScale.bandwidth())
        });

        canvasSelection.on("mousemove", () => {
            const mouse = d3_mouse(canvas);
            const mouseX = mouse[0];
            const mouseY = mouse[1];

            const y = yScale.invert(mouseY);

            const mouseViewportX = d3_event.clientX;
            const mouseViewportY = d3_event.clientY;
            
            PubSub.publish("tooltip", {
                x: mouseViewportX,
                y: mouseViewportY,
                content: "Test " + y
            });
        });

        canvasSelection.on("mouseout", () => {
            PubSub.publish("tooltip", {
                x: null,
                y: null,
                content: null
            });
        })

        return (() => teardownCanvas(canvasRef));
    });

    console.log("TrackRowInfo.render");
    return (
        <div>
            <canvas
                ref={canvasRef}
                style={{
                    position: 'absolute',
                    top: `${top}px`,
                    left: `${left}px`, 
                    width: `${width}px`,
                    height: `${height}px`,
                }}
            />
        </div>
    );
}