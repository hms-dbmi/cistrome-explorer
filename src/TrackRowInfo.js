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
 * @prop {array} infoAttributes Array of attribute names.
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 */
export default function TrackRowInfo(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight, 
        rowInfo, 
        infoAttributes,
        rowInfoPosition
    } = props;

    // Dimensions
    const top = trackY;
    const colWidth = 15;
    const xMargin = 60;
    const xMarginInitial = 5;
    const width = xMargin * infoAttributes.length;
    const height = trackHeight;

    // Scales
    const xScale = d3_scaleThreshold();
    const yScale = vega_scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);
    
    // Stores valueScale, left position, attribute name for each attribute
    let vizRecipes = [];

    // Viz recipes condition on left vs. right positioning:
    let left, xScaleDomain = [], xScaleRange = [];
    if(rowInfoPosition === "left") {
        left = trackX - xMarginInitial - width;

        xScaleRange.push(`margin-${infoAttributes.length}`);
        for(let i = infoAttributes.length - 1; i >= 0; i--){    // First attribute is shown on the 'right-most' side
            const attribute = infoAttributes[i];

            const dimLeft = (colWidth + xMargin) * (infoAttributes.length - 1 - i);
            const colorScale = d3_scaleOrdinal()
                .domain(Array.from(new Set(rowInfo.map(d => d[attribute]))))
                .range(d3_schemeCategory10);
            
            vizRecipes.push({left: dimLeft, colorScale, attribute});

            // [secondaryLeft, secondaryLeft+colWidth, primaryLeft, primaryLeft+colWidth, ... width]
            xScaleDomain.push(dimLeft, dimLeft + colWidth);
            // ["m-n", ... "d-1", "m-1", "d-0", "m-0"]
            xScaleRange.push(`dimension-${i}`, `margin-${i}`);
        }
        xScaleDomain.push(width);

    } else if(rowInfoPosition === "right") {
        left = trackWidth + trackX + xMarginInitial;

        for(let i = 0; i < infoAttributes.length; i++){
            const attribute = infoAttributes[i];

            const dimLeft = (colWidth + xMargin) * i;
            const colorScale = d3_scaleOrdinal()
                .domain(Array.from(new Set(rowInfo.map(d => d[attribute]))))
                .range(d3_schemeCategory10);

            vizRecipes.push({left: dimLeft, colorScale, attribute});

            // [primaryLeft, primaryLeft+colWidth, secondaryLeft, secondaryLeft+colWidth, ... width]
            xScaleDomain.push(dimLeft, dimLeft + colWidth);
            // ["m-0", "d-0", "m-1", "d-1", ..., "m-n"]
            xScaleRange.push(`margin-${i}`, `dimension-${i}`);
        }
        xScaleDomain.push(width);
        xScaleRange.push(`margin-${infoAttributes.length}`);
    }
    
    // X-axis scale for mouse events
    xScale
        .domain(xScaleDomain)
        .range(xScaleRange);
    
    // Canvas
    const canvasRef = useRef();
    useEffect(() => {
        const { canvas, context, canvasSelection } = setupCanvas(canvasRef);
        context.clearRect(0, 0, width, height);

        // Draw metadata value colors.
        vizRecipes.forEach(recipe => {
            const {left, colorScale, attribute} = recipe;
            rowInfo.forEach((d, i) => {
                context.fillStyle = colorScale(d[attribute]);
                context.fillRect(left, yScale(i), colWidth, yScale.bandwidth());
            });
        });

        canvasSelection.on("mousemove", () => {
            const mouse = d3_mouse(canvas);
            const mouseX = mouse[0];
            const mouseY = mouse[1];

            const y = yScale.invert(mouseY);
            const x = xScale(mouseX);

            let xVal;
            if(x.includes("dimension")){
                const dimIndex = parseInt(x.split("-")[1], 10);
                xVal = rowInfo[y][infoAttributes[dimIndex]];
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