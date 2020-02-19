import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { drawVisTitle } from "./utils/vis.js";

import TrackRowInfoControl from './TrackRowInfoControl.js';

export const margin = 5;

/**
 * Component for visualization of row info quantitative or nominal attribute values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {string} sortOrder The order of sort applied in this track. If sort is not applied, then null.
 * @prop {string} viewId The viewId for the horizontal-multivec track.
 * @prop {string} trackId The trackId for the horizontal-multivec track.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisBar(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        viewId,
        trackId,
        rowInfo,
        sortOrder,
        onSortRows,
        onSearchRows,
        drawRegister,
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();
    const [mouseX, setMouseX] = useState(null);

    // Data, layouts and styles
    const { field, type } = fieldInfo;
    const isNominal = type === "nominal";

    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);
    const rowHeight = yScale.bandwidth();

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        drawVisTitle(field, { two, isLeft, isNominal, width, sortOrder });

        const textAreaWidth = isNominal ? width - 20 : 20;
        const barAreaWidth = width - textAreaWidth;
        const minTrackWidth = 40;
        const isTextLabel = width > minTrackWidth;
        const fontSize = 10;
        
        // Scales
        const valueExtent = [0, d3.extent(rowInfo.map(d => d[field]))[1]];   // Zero baseline
       
        const xScale = d3.scaleLinear()
            .domain(valueExtent)
            .range([0, barAreaWidth]);

        const colorScale = isNominal ? 
            d3.scaleOrdinal()
                .domain(Array.from(new Set(rowInfo.map(d => d[field]))).sort()) 
                .range(d3.schemeSet3) : 
            d3.scaleLinear()
                .domain(valueExtent)
                .range([0, 1]);
    

        // Render visual components for each row (i.e., bars and texts).
        const textAlign = isLeft ? "end" : "start";
        let aggregateStartIdx = -1, sameCategoriesNearby = 1;

        rowInfo.forEach((d, i) => {
            // To aggregate bars, check if there is a same category on the next row.
            if(type === "nominal" && i + 1 < rowInfo.length && d[field] === rowInfo[i+1][field]) {
                if(aggregateStartIdx === -1) {
                    aggregateStartIdx = i;
                }
                sameCategoriesNearby++;
                return;
            }

            const barTop = aggregateStartIdx !== -1 ? yScale(aggregateStartIdx) : yScale(i);
            const barHeight = rowHeight * sameCategoriesNearby;
            const barWidth = isNominal ? barAreaWidth : xScale(d[field]);        
            const barLeft = (isLeft ? width - barWidth : 0);
            const textLeft = (isLeft ? width - barWidth - margin : barWidth + margin);
            const color = isNominal ? colorScale(d[field]) : 
                d3.interpolateViridis(colorScale(d[field]));

            const rect = two.makeRect(barLeft, barTop, barWidth, barHeight);
            rect.fill = color;

            // Render text labels when the space is enough.
            if(barHeight >= fontSize && isTextLabel){
                const text = two.makeText(textLeft, barTop + barHeight/2, textAreaWidth, barHeight, d[field]);
                text.fill = d3.hsl(color).darker(3);
                text.fontsize = fontSize;
                text.align = textAlign;
                text.baseline = "middle";
                text.overflow = "ellipsis";
            }

            aggregateStartIdx = -1;
            sameCategoriesNearby = 1;
        });

        two.update();
        return two.teardown;
    });
    
    drawRegister("TrackRowInfoVisBar", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            let fieldVal;
            if(y !== undefined){
                setMouseX(true);
                fieldVal = rowInfo[y][field];
            } else {
                setMouseX(null);
                destroyTooltip();
                return;
            }

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;
            
            PubSub.publish(EVENT.TOOLTIP, {
                x: mouseViewportX,
                y: mouseViewportY,
                content: `${field}: ${fieldVal}`
            });
        });

        // Handle mouse leave.
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => setMouseX(null));

        // Clean up.
        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, rowInfo]);

    return (
        <div
            ref={divRef}
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
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    position: 'relative'
                }}
            />
            <TrackRowInfoControl
                isLeft={isLeft}
                isVisible={mouseX !== null}
                fieldInfo={fieldInfo}
                searchTop={top}
                searchLeft={left}
                onSortRows={onSortRows}
                onSearchRows={onSearchRows}
            />
        </div>
    );
}