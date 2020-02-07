import React, { useEffect, useRef, useCallback } from 'react';
import range from 'lodash/range';
import d3 from './utils/d3.js';
import Two from './utils/two.js';

import { EVENT } from './constants.js';

import './TrackRowInfo.scss';
import { visualizationTrack } from './visualizationTrack.js';

function destroyTooltip() {
    PubSub.publish(EVENT.TOOLTIP, {
        x: null,
        y: null,
        content: null
    });
}

/**
 * Component for visualization of row info attribute values.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {array} infoAttributes Array of JSON object, one object for the names and types of each attribute.
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 * @prop {function} register The function for child components to call to register their draw functions.
 */
export default function TrackRowInfo(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight, 
        rowInfo, 
        infoAttributes,
        rowInfoPosition,
        register
    } = props;

    // Dimensions
    const isLeft = rowInfoPosition === "left";
    const top = trackY;
    const unitWidth = 100;
    const width = unitWidth * infoAttributes.length;
    const margin = 5;
    const height = trackHeight;
    const left = isLeft ? trackX - margin - width : trackX + trackWidth + margin;

    // Render canvas
    const canvasRef = useRef();

    // Determine position of each dimension.
    let xDomain = [], xRange = [];
    for(let i = 0; i < infoAttributes.length; i++) {
        const attribute = isLeft ? infoAttributes[infoAttributes.length - i - 1] : infoAttributes[i];
        let currentLeft = unitWidth * i;

        // Domain and range for mouse event
        xDomain.push(currentLeft + unitWidth);
        xRange.push(attribute.name);
    }
    
    // Scales
    const xScale = d3.scaleThreshold()
        .domain(xDomain)
        .range(xRange);
    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });
      
        for(let i = 0; i < infoAttributes.length; i++) {
            const attribute = isLeft ? infoAttributes[infoAttributes.length - i - 1] : infoAttributes[i];
            let currentLeft = unitWidth * i;

            visualizationTrack({
                two, 
                left: currentLeft, top: 0, width: unitWidth, height: height,
                rowInfo,
                attribute,
                isLeft
            });
        }

        two.update();
        return two.teardown;
    }, [width, height]);

    register("TrackRowInfo", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
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
                content: `${x}: ${xVal}`
            });
        });

        d3.select(canvas).on("mouseout", destroyTooltip)

        return teardown;
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