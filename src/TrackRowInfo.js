import React, { useState, useEffect, useRef, useCallback } from 'react';
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
        rowInfoAttributes,
        rowInfoPosition,
        register
    } = props;

    const [mouseX, setMouseX] = useState(-1);

    // Dimensions
    const isLeft = rowInfoPosition === "left";
    const top = trackY;
    const unitWidth = 100;
    const width = unitWidth * rowInfoAttributes.length;
    const height = trackHeight;
    const left = isLeft ? trackX - width : trackX + trackWidth;

    // Render canvas
    const canvasRef = useRef();

    // Determine position of each dimension.
    let trackProps = [], xDomain = [], xRange = [];
    rowInfoAttributes.forEach((attribute, i) => {
        const fieldInfo = isLeft ? rowInfoAttributes[rowInfoAttributes.length - i - 1] : attribute;
        let currentLeft = unitWidth * i;

        trackProps.push({
            left: currentLeft, top: 0, width: unitWidth, height: height,
            rowInfo,
            fieldInfo,
            isLeft
        });

        // Domain and range for mouse event.
        xDomain.push(currentLeft + unitWidth);
        xRange.push(fieldInfo.field);
    });
    
    // Scales
    const xScale = d3.scaleThreshold()
        .domain(xDomain)
        .range(xRange);
    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);

    // Render each track.
    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });
    
        trackProps.forEach(d => visualizationTrack({...d, two}));

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
            if(y !== undefined && x !== undefined){
                setMouseX(rowInfoAttributes.map(d => d.field).indexOf(x));
                xVal = rowInfo[y][x];
            } else {
                setMouseX(-1);
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
        d3.select(canvas).on("mouseout", () => destroyTooltip());
        return teardown;
    });

    // Make small control panels for each track.
    let trackControls = trackProps.map((d, i) => (
        <div 
            key={i}
            className={"cistrome-track-control"}
            style={{
                top: `${d.top + 2}px`,
                left: `${d.left + 2}px`, 
                width: `${40}px`,
                height: `${20}px`,
                visibility: mouseX === i ? "visible" : "hidden"
            }}
        >
            <svg
                className="cistrome-track-control-button-left"
                onClick={() => console.log("button clicked!")}
            >
                <title>Sort rows in ascending order</title>
                <use xlinkHref="#chevron_up" />
            </svg>
            <svg
                className="cistrome-track-control-button-right"
                onClick={() => console.log("button clicked!")}
            >
                <title>Sort rows in descending order</title>
                <use xlinkHref="#chevron_down" />
            </svg>
        </div>
    ), this);

    function onMouseLeave() {
        setMouseX(-1);
        destroyTooltip();
    };

    return (
        <div className="cistrome-hgw-child"
            onMouseLeave={onMouseLeave}
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
            {trackControls}
        </div>
    );
}