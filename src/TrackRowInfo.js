import React, { useEffect, useRef, useCallback } from 'react';
import range from 'lodash/range';
import d3 from './utils/d3.js';
import Two from './utils/two.js';
import { teardownCanvas } from './utils/canvas.js';

import { EVENT } from './constants.js';

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
        rowInfoPosition,
        register
    } = props;

    // Dimensions
    const top = trackY;
    const colWidth = 10;    // width of stacked bars
    const xMargin = 60;     // width of text area
    const xMarginInitial = 5;
    const xGap = 5; // gap between bars and text
    const width = (colWidth + xMargin) * infoAttributes.length;
    const height = trackHeight;
    const titleFontSize = 12;
    const fontSize = 10;

    // Scales
    const xScale = d3.scaleThreshold();
    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);
    
    // Stores recipes to visualize each attribute using texts and color bars
    let vizRecipes = [];

    // Viz recipes depending on left or right positioning
    let left, xScaleDomain = [], xScaleRange = [];
    if(rowInfoPosition === "left") {
        left = trackX - xMarginInitial - width;

        xScaleRange.push(`margin-${infoAttributes.length}`);
        for(let i = infoAttributes.length - 1; i >= 0; i--) {    // First attribute is shown on the 'right-most' side
            const attribute = infoAttributes[i];

            const dimLeft = xMargin + (colWidth + xMargin) * (infoAttributes.length - 1 - i);
            const colorScale = d3.scaleOrdinal()
                .domain(Array.from(new Set(rowInfo.map(d => d[attribute]))))
                .range(d3.schemeSet3);
            
            vizRecipes.push({
                titleLeft: dimLeft - xMargin + xGap, 
                titleRotate: -Math.PI/2,
                titleTextAlign: "right", 
                barLeft: dimLeft, 
                labelLeft: dimLeft - xGap, 
                textAlign: "right", 
                colorScale, 
                attribute
            });

            // [secondaryLeft, secondaryLeft+colWidth, primaryLeft, primaryLeft+colWidth, ... width]
            xScaleDomain.push(dimLeft, dimLeft + colWidth);
            // ["m-n", ... "d-1", "m-1", "d-0", "m-0"]
            xScaleRange.push(`dimension-${i}`, `margin-${i}`);
        }
        xScaleDomain.push(width);

    } else if(rowInfoPosition === "right") {
        left = trackWidth + trackX + xMarginInitial;

        for(let i = 0; i < infoAttributes.length; i++) {
            const attribute = infoAttributes[i];

            const dimLeft = (colWidth + xMargin) * i;
            const colorScale = d3.scaleOrdinal()
                .domain(Array.from(new Set(rowInfo.map(d => d[attribute]))))
                .range(d3.schemeSet3);

            vizRecipes.push({
                titleLeft: dimLeft + colWidth + xMargin - xGap, 
                titleRotate: Math.PI/2,
                titleTextAlign: "left", 
                barLeft: dimLeft, 
                labelLeft: dimLeft + colWidth + xGap,
                textAlign: "left", 
                colorScale, 
                attribute
            });

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

    const rowHeight = yScale.bandwidth();
    
    // Canvas
    const canvasRef = useRef();

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });
        
        // Show metadata values with visual elements
        vizRecipes.forEach(recipe => {
            const {
                titleLeft, titleRotate, titleTextAlign,
                barLeft, labelLeft, textAlign, colorScale, 
                attribute
            } = recipe;

            // Draw a title of each dimension
            const title = two.makeText(titleLeft, 0, rowHeight, colWidth, `attribute: ${attribute}`);
            title.fill = "#9A9A9A";
            title.fontsize = titleFontSize;
            title.align = titleTextAlign;
            title.baseline = "top";
            title.rotation = titleRotate;

            // Draw color bars and text labels
            rowInfo.forEach((d, i) => {
                const rect = two.makeRect(barLeft + colWidth/2, yScale(i) + rowHeight/2, colWidth, rowHeight);
                rect.fill = colorScale(d[attribute]);

                if(rowHeight >= fontSize){
                    const text = two.makeText(labelLeft, yScale(i) + rowHeight/2, colWidth, rowHeight, d[attribute])
                    text.fill = d3.hsl(colorScale(d[attribute])).darker(3);
                    text.fontsize = fontSize;
                    text.align = textAlign;
                    text.baseline = "middle";
                }
            });
        });

        two.update();
    }, [width, height, vizRecipes]);

    register("TrackRowInfo", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        
        draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const mouse = d3.mouse(canvas);
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

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;
            
            PubSub.publish(EVENT.TOOLTIP, {
                x: mouseViewportX,
                y: mouseViewportY,
                content: `Value: ${xVal}`
            });
        });

        d3.select(canvas).on("mouseout", destroyTooltip)

        return (() => teardownCanvas(canvas));
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