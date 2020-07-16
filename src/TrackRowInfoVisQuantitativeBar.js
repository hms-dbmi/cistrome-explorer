import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { TooltipContent, destroyTooltip } from "./Tooltip.js";
import { drawVisTitle } from "./utils/vis.js";

import TrackRowInfoControl from './TrackRowInfoControl.js';
import { rgbToHex, generateNextUniqueColor } from "./utils/color.js";
import { getRetinaRatio } from './utils/canvas.js';
import { modifyItemInArray } from "./utils/array.js";
import { getAggregatedValue } from "./utils/aggregate.js";

export const margin = 5;

/**
 * Component for visualization of row info quantitative or nominal attribute values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {boolean} isShowControlButtons Determine if control buttons should be shown.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {string} titleSuffix The suffix of a title, information about sorting and filtering status.
 * @prop {object} sortInfo The options for sorting rows of the field used in this track.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a highlight interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisQuantitativeBar(props) {
    const {
        left, top, width, height,
        field, type, title, aggFunction,
        isLeft,
        isShowControlButtons,
        rowInfo,
        transformedRowInfo,
        titleSuffix,
        sortInfo,
        filterInfo,
        onAddTrack,
        onSortRows,
        onHighlightRows,
        onFilterRows,
        drawRegister,
    } = props;

    const divRef = useRef();
    const axisRef = useRef();
    const canvasRef = useRef();
    const hiddenCanvasRef = useRef();

    // Data, layouts and styles
    const isStackedBar = Array.isArray(field);
    const axisHeight = 30;
    const textAreaWidth = 20;
    const barAreaWidth = width - textAreaWidth;
    const minTrackWidth = 40;
    const fontSize = 10;
    const aggValue = (d, f) => getAggregatedValue(d, f, "quantitative", aggFunction);
    const numberFormatShort = d3.format(".0f");
    const numberFormatLong = d3.format(".2f");

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
                value: aggValue(d, field),
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
        const isTextLabel = width > minTrackWidth;

        if(isStackedBar) {
            // Scales
            const valueExtent = [0, d3.extent(transformedRowInfo.map(d => {
                let sum = 0;
                field.forEach(f => sum += aggValue(d, f));
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
                    const barWidth = xScale(aggValue(d, f));
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
                    field.forEach(f => sum += aggValue(d, f));
                    let textLeft = (isLeft ? width - xScale(sum) - margin : xScale(sum) + margin);
                    const text = two.makeText(textLeft, barTop + rowHeight/2, textAreaWidth, rowHeight, numberFormatShort(sum));
                    text.fill = "black";
                    text.fontsize = fontSize;
                    text.align = textAlign;
                    text.baseline = "middle";
                    text.overflow = "clip";
                }
            });

        } else {
            // Scales
            const valueExtent = [0, d3.extent(transformedRowInfo.map(d => aggValue(d, field)))[1]];   // Zero baseline
            xScale = d3.scaleLinear()
                .domain(valueExtent)
                .range([0, barAreaWidth]);
            const colorScale = d3.scaleLinear()
                    .domain(valueExtent)
                    .range([0, 1]);

            // Render visual components for each row (i.e., bars and texts).
            const textAlign = isLeft ? "end" : "start";
            transformedRowInfo.forEach((d, i) => {
                const value = numberFormatShort(aggValue(d, field));
                const barTop = yScale(i);
                const barWidth = xScale(value);
                const barLeft = (isLeft ? width - barWidth : 0);
                const textLeft = (isLeft ? width - barWidth - margin : barWidth + margin);
                const infoForMouseEvent = colorToInfo.find(d => d.field === field && d.rowIndex === i);
                const color = isHidden ? infoForMouseEvent.uniqueColor : d3.interpolateViridis(colorScale(value));;

                const rect = two.makeRect(barLeft, barTop, barWidth, rowHeight);
                rect.fill = color;

                // Render text labels when the space is enough.
                if(rowHeight >= fontSize && isTextLabel) {
                    const text = two.makeText(textLeft, barTop + rowHeight/2, textAreaWidth, rowHeight, numberFormatShort(value));
                    text.fill = d3.hsl(color).darker(3);
                    text.fontsize = fontSize;
                    text.align = textAlign;
                    text.baseline = "middle";
                    text.overflow = "clip";
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

        drawVisTitle(title, { two, isLeft, width, height, titleSuffix });

        two.update();
        return two.teardown;
    });

    const drawAxis = useCallback((domElement) => {
        d3.select(domElement).selectAll("*").remove();

        const axisScale = isLeft ? xScale.domain(xScale.domain().reverse()) : xScale;
        const axis = d3.axisBottom(axisScale)
            .ticks(Math.ceil(barAreaWidth / 40));
        
        d3.select(domElement)
            .attr("width", width)
            .attr("height", axisHeight)
            .append("g")
                .attr("transform", `translate(${isLeft ? textAreaWidth - 1 : 1}, 0)`)
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
                        title={(isStackedBar ? `${title}: ${hoveredInfo.field}` : title)}
                        value={numberFormatLong(hoveredInfo.value)}
                        color={hoveredInfo.color}
                    />
                });
            } else {
                destroyTooltip();
            }            
        });

        // Handle mouse enter and leave.
        d3.select(canvas).on("mouseout", destroyTooltip);

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
            style={{
                top: `${top}px`,
                position: 'relative',
                width: `${width}px`,
                height: `${height}px`,
            }}
        >
            <canvas
                ref={canvasRef}
                style={{
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    position: 'absolute'
                }}
            />
            <canvas
                ref={hiddenCanvasRef}
                className="hm-hidden"
                style={{
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    position: 'absolute'
                }}
            />
            <svg ref={axisRef} 
                style={{ 
                    pointerEvents: "none",
                    position: "absolute"
                }}
            />
            <TrackRowInfoControl
                isLeft={isLeft}
                isVisible={isShowControlButtons}
                field={field}
                type={type}
                title={title}
                aggFunction={aggFunction}
                searchTop={top + axisHeight}
                searchLeft={left}
                sortAsceButtonHighlit={sortInfo && sortInfo.order === "ascending"}
                sortDescButtonHighlit={sortInfo && sortInfo.order === "descending"}
                filterButtonHighlit={filterInfo !== undefined}
                onSortRows={onSortRows}
                onHighlightRows={onHighlightRows}
                onFilterRows={onFilterRows}
                filterInfo={filterInfo}
                rowInfo={rowInfo}
                transformedRowInfo={transformedRowInfo}
            />
        </div>
    );
}