import React, { useEffect, useRef } from 'react';
import range from 'lodash/range';
import d3 from './d3.js';

import { EVENT } from './constants.js';
import { setupCanvas, teardownCanvas } from './utils-canvas.js';
import { verticalBarTrack } from './VerticalBarTrack.js';

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
 * @prop {array} infoAttributes Array of JSON object, one object for the names and types of each attribute.
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
    const isLeft = rowInfoPosition === "left";
    const top = trackY;
    const barWidth = 10;
    const textWidth = 60;
    const margin = 5;
    const width = (barWidth + textWidth) * infoAttributes.length;
    const height = trackHeight;
    const left = isLeft ? trackX - margin - width : trackX + trackWidth + margin;

    // Render canvas
    const canvasRef = useRef();
    useEffect(() => {
        const { canvas, context, canvasSelection } = setupCanvas(canvasRef);
        context.clearRect(0, 0, width, height);

        // Determin position of each dimension and render it
        let xDomain = [], xRange = [];
        for(let i = 0; i < infoAttributes.length; i++) {
            const attribute = isLeft ? infoAttributes[infoAttributes.length - i - 1] : infoAttributes[i];
            let currentLeft = (barWidth + textWidth) * i;

            verticalBarTrack({
                ref: context, 
                left: currentLeft, top: 0, width: barWidth + textWidth, height: height,
                rowInfo,
                attribute,
                isLeft, isCanvas: true
            });

            // Domain and range for mouse event
            if(isLeft) {
                xDomain.push(currentLeft + textWidth, currentLeft + textWidth + barWidth);
                xRange.push("m", attribute.name);
            } else {
                xDomain.push(currentLeft, currentLeft + barWidth);
                xRange.push("m", attribute.name);
            }
        }
        xDomain.push(width);
        xRange.push("m");
        
        // Scales
        const xScale = d3.scaleThreshold()
            .domain(xDomain)
            .range(xRange);
        const yScale = d3.scaleBand()
            .domain(range(rowInfo.length))
            .range([0, height]);

        canvasSelection.on("mousemove", () => {
            const mouse = d3.mouse(canvas);
            const mouseX = mouse[0];
            const mouseY = mouse[1];

            const y = yScale.invert(mouseY);
            const x = xScale(mouseX);
            let xVal;
            if(x !== "m"){
                xVal = rowInfo[y][x];
            } else {
                destroyTooltip();
                return;
            }

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;
            
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