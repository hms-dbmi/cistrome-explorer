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
 * @prop {object[]} transformedRowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {string} titleSuffix The suffix of a title, information about sorting and filtering status.
 * @prop {string} viewId The viewId for the horizontal-multivec track.
 * @prop {string} trackId The trackId for the horizontal-multivec track.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisQuantitativeBar(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        viewId,
        trackId,
        transformedRowInfo,
        titleSuffix,
        onSortRows,
        onSearchRows,
        onFilter,
        drawRegister,
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();
    const [mouseX, setMouseX] = useState(null);

    // Data, layouts and styles
    const { field } = fieldInfo;
    const isStackedBar = Array.isArray(field);

    const yScale = d3.scaleBand()
        .domain(range(transformedRowInfo.length))
        .range([0, height]);
    const rowHeight = yScale.bandwidth();

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        const titleText = isStackedBar ? field.join(" + ") : field;
        drawVisTitle(titleText, { two, isLeft, isNominal: false, width, titleSuffix });

        const textAreaWidth = 20;
        const barAreaWidth = width - textAreaWidth;
        const minTrackWidth = 40;
        const isTextLabel = width > minTrackWidth;
        const fontSize = 10;
       
        if(isStackedBar) {
            // Scales
            const valueExtent = [0, d3.extent(transformedRowInfo.map(d => {
                let sum = 0;
                field.forEach(f => sum += d[f]);
                return sum;
            }))[1]];   // Zero baseline
            const xScale = d3.scaleLinear()
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

                field.forEach((f, i) => {
                    const barWidth = xScale(d[f]);    
                    const color = colorScale(f);

                    currentBarLeft += (isLeft ? -barWidth : 0);

                    const rect = two.makeRect(currentBarLeft, barTop, barWidth, rowHeight);
                    rect.fill = color;

                    currentBarLeft += (isLeft ? 0 : barWidth);
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
            const xScale = d3.scaleLinear()
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
                const color = d3.interpolateViridis(colorScale(d[field]));

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
            });
        }

        two.update();
        return two.teardown;
    });
    
    drawRegister("TrackRowInfoVisQuantitativeBar", draw);

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
                if(Array.isArray(field)){
                    fieldVal = 0;
                    field.forEach(f => fieldVal += transformedRowInfo[y][f]);

                } else {
                    fieldVal = transformedRowInfo[y][field];
                }
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
                content: Array.isArray(field) ? `${field.join(' + ')}: ${fieldVal}` : `${field}: ${fieldVal}`
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
    }, [top, left, width, height, transformedRowInfo]);

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
                onFilter={onFilter}
            />
        </div>
    );
}