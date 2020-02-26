import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { TooltipContent } from "./Tooltip.js";
import { drawVisTitle } from "./utils/vis.js";

import TrackRowInfoControl from './TrackRowInfoControl.js';
import { rgbToHex, generateNextUniqueColor } from "./utils/color.js";
import { getRetinaRatio } from './utils/canvas.js';
import { modifyItemInArray } from "./utils/array.js";

export const margin = 5;

/**
 * Component for visualization of row info quantitative or nominal attribute values.
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
export default function TrackRowInfoVisQuantitativeBar(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        transformedRowInfo,
        titleSuffix,
        onSortRows,
        onSearchRows,
        onFilterRows,
        drawRegister,
    } = props;

    const divRef = useRef();
    const axisRef = useRef();
    const canvasRef = useRef();
    const hiddenCanvasRef = useRef();
    const [isMouseHover, setIsMouseHover] = useState(null);

    // Data, layouts and styles
    const { field } = fieldInfo;
    const isStackedBar = Array.isArray(field);
    const axisHeight = 30;
    const textAreaWidth = 20;
    const barAreaWidth = width - textAreaWidth;
    const minTrackWidth = 40;
    const fontSize = 10;

    let xScale = d3.scaleLinear();
    const yScale = d3.scaleBand()
        .domain(range(transformedRowInfo.length))
        .range([0, height]);
    const rowHeight = yScale.bandwidth();

    // Array to store information for mouse events, such as unique color.
    let colorToInfo = [];
    let cnt = 1, fields = isStackedBar ? field : [field];
    
    fields.forEach(field => {
        transformedRowInfo.forEach((d, i) => {
            const uniqueColor = generateNextUniqueColor(cnt++);
            colorToInfo.push({
                uniqueColor,
                field,
                value: d[field],
                rowIndex: i,
                color: null // This property is determined when first rendered.
            });
            cnt += 1;
        });
    });

    const draw = useCallback((domElement, isHidden) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        const titleText = isStackedBar ? field.join(" + ") : field;

        const isTextLabel = width > minTrackWidth;
       
        if(isStackedBar) {
            // Scales
            const valueExtent = [0, d3.extent(transformedRowInfo.map(d => {
                let sum = 0;
                field.forEach(f => sum += d[f]);
                return sum;
            }))[1]];   // Zero baseline
            xScale = xScale
                .domain(valueExtent)
                .range([0, barAreaWidth]);
            const colorScale = d3.scaleOrdinal()
                .domain(Array.from(new Set(field)).sort())
                .range(d3.schemeTableau10);

            // Render visual components for each row (i.e., bars and texts).
            const textAlign = isLeft ? "end" : "start";
            transformedRowInfo.forEach((d, i) => {
                const barTop = yScale(i);
                let currentBarLeft = (isLeft ? width : 0);

                field.forEach(f => {
                    const barWidth = xScale(d[f]);
                    const infoForMouseEvent = colorToInfo.find(d => d.field === f && d.rowIndex === i);
                    const color = isHidden ? infoForMouseEvent.uniqueColor : colorScale(f);
                    
                    currentBarLeft += (isLeft ? -barWidth : 0);

                    const rect = two.makeRect(currentBarLeft, barTop, barWidth, rowHeight);
                    rect.fill = color;

                    currentBarLeft += (isLeft ? 0 : barWidth);

                    // Add other vis properties to colorToInfo for mouse events.
                    if(!isHidden) {
                        colorToInfo = modifyItemInArray(colorToInfo, colorToInfo.indexOf(infoForMouseEvent), {
                            ...infoForMouseEvent, 
                            color
                        });
                    }
                });

                // Render text labels when the space is enough.
                if(rowHeight >= fontSize && isTextLabel) {
                    let sum = 0;
                    field.forEach(f => sum += d[f]);
                    let textLeft = (isLeft ? width - xScale(sum) - margin : xScale(sum) + margin);
                    const text = two.makeText(textLeft, barTop + rowHeight/2, textAreaWidth, rowHeight, sum);
                    text.fill = "black";
                    text.fontsize = fontSize;
                    text.align = textAlign;
                    text.baseline = "middle";
                    text.overflow = "ellipsis";
                }
            });
        } else {
            // Scales
            const valueExtent = [0, d3.extent(transformedRowInfo.map(d => d[field]))[1]];   // Zero baseline
            xScale = d3.scaleLinear()
                .domain(valueExtent)
                .range([0, barAreaWidth]);
            const colorScale = d3.scaleLinear()
                    .domain(valueExtent)
                    .range([0, 1]);

            // Render visual components for each row (i.e., bars and texts).
            const textAlign = isLeft ? "end" : "start";
            transformedRowInfo.forEach((d, i) => {
                const barTop = yScale(i);
                const barWidth = xScale(d[field]);        
                const barLeft = (isLeft ? width - barWidth : 0);
                const textLeft = (isLeft ? width - barWidth - margin : barWidth + margin);
                const infoForMouseEvent = colorToInfo.find(d => d.field === field && d.rowIndex === i);
                const color = isHidden ? infoForMouseEvent.uniqueColor : d3.interpolateViridis(colorScale(d[field]));;

                const rect = two.makeRect(barLeft, barTop, barWidth, rowHeight);
                rect.fill = color;

                // Render text labels when the space is enough.
                if(rowHeight >= fontSize && isTextLabel) {
                    const text = two.makeText(textLeft, barTop + rowHeight/2, textAreaWidth, rowHeight, d[field]);
                    text.fill = d3.hsl(color).darker(3);
                    text.fontsize = fontSize;
                    text.align = textAlign;
                    text.baseline = "middle";
                    text.overflow = "ellipsis";
                }

                // Add other vis properties to colorToInfo for mouse events.
                if(!isHidden) {
                    colorToInfo = modifyItemInArray(colorToInfo, colorToInfo.indexOf(infoForMouseEvent), {
                        ...infoForMouseEvent, 
                        color
                    });
                }
            });
        }

        drawVisTitle(titleText, { two, isLeft, width, height, titleSuffix });

        two.update();
        return two.teardown;
    });

    const drawAxis = useCallback((domElement) => {
        d3.select(domElement).selectAll("*").remove();

        const axisScale = isLeft ? xScale.domain(xScale.domain().reverse()) : xScale;
        const axis = d3.axisTop(axisScale)
            .ticks(Math.ceil(barAreaWidth / 40));
        
        d3.select(domElement)
            .attr("width", width)
            .attr("height", axisHeight)
            .append("g")
                .attr("transform", `translate(${isLeft ? textAreaWidth - 1 : 1}, ${axisHeight - margin})`)
                .call(axis);
        
        d3.select(domElement)
            .selectAll("text")
                .attr("transform", `translate(${isLeft ? -3 : 3}, 0)`);
        
        return () => { /* Teardown */ };
    });
    
    drawRegister("TrackRowInfoVisQuantitativeBar", draw);
    drawRegister("TrackRowInfoVisQuantitativeBarAxis", drawAxis);

    useEffect(() => {
        const canvas = canvasRef.current;
        const hiddenCanvas = hiddenCanvasRef.current;
        const svg = axisRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);
        const teardownHidden = draw(hiddenCanvas, true);
        const teardownSvg = drawAxis(svg);

        d3.select(canvas).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);
            
            const hiddenContext = hiddenCanvasRef.current.getContext('2d');
            const ratio = getRetinaRatio(hiddenContext);
            const uniqueColor = rgbToHex(hiddenContext.getImageData(mouseX * ratio, mouseY * ratio, 1, 1).data);
            const hoveredInfo = colorToInfo.find(d => d.uniqueColor === uniqueColor);

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;
            
            if(hoveredInfo) {
                PubSub.publish(EVENT.TOOLTIP, {
                    x: mouseViewportX,
                    y: mouseViewportY,
                    content: <TooltipContent 
                        title={hoveredInfo.field}
                        value={hoveredInfo.value}
                        color={hoveredInfo.color}
                    />
                });
            } else {
                destroyTooltip();
            }            
        });

        // Handle mouse enter and leave.
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseenter", () => setIsMouseHover(true));
        d3.select(div).on("mouseleave", () => setIsMouseHover(null));

        // Clean up.
        return () => {
            teardown();
            teardownHidden();
            teardownSvg();
            d3.select(div).on("mouseenter", null);
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, transformedRowInfo]);

    return (
        <div
            ref={divRef}
            className="cistrome-hgw-child"
            style={{
                top: `${top - axisHeight}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height + axisHeight}px`,
            }}
        >
            <canvas
                ref={canvasRef}
                style={{
                    top: axisHeight,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    position: 'absolute'
                }}
            />
            <canvas
                ref={hiddenCanvasRef}
                className="chw-hidden-canvas"
                style={{
                    top: axisHeight,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    position: 'absolute'
                }}
            />
            <svg ref={axisRef}
                style={{
                    color: "#9A9A9A"
                }}
            />
            <TrackRowInfoControl
                isLeft={isLeft}
                isVisible={isMouseHover !== null}
                fieldInfo={fieldInfo}
                searchTop={top + axisHeight}
                searchLeft={left}
                onSortRows={onSortRows}
                onSearchRows={onSearchRows}
                onFilterRows={onFilterRows}
            />
        </div>
    );
}