import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT, CONTEXT_MENU_TYPE } from "./utils/constants.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { drawVisTitle } from "./utils/vis.js";

import TrackRowInfoControl from './TrackRowInfoControl.js';
import { TooltipContent } from "./Tooltip.js";

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
export default function TrackRowInfoVisNominalBar(props) {
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
    const canvasRef = useRef();
    const [mouseX, setMouseX] = useState(null);
    const [hoverValue, setHoverValue] = useState();

    // Data, layouts and styles
    const { field } = fieldInfo;

    const yScale = d3.scaleBand()
        .domain(range(transformedRowInfo.length))
        .range([0, height]);
    const rowHeight = yScale.bandwidth();
    const colorScale = d3.scaleOrdinal()
        .domain(Array.from(new Set(transformedRowInfo.map(d => d[field]))).sort()) 
        .range(d3.schemeTableau10);

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        const titleText = Array.isArray(field) ? field.join(" + ") : field;
        
        const textAreaWidth = width - 20;
        const barAreaWidth = width - textAreaWidth;
        const minTrackWidth = 40;
        const isTextLabel = width > minTrackWidth;
        const fontSize = 10;

        // Render visual components for each row (i.e., bars and texts).
        const textAlign = isLeft ? "end" : "start";
        let aggregateStartIdx = -1, sameCategoriesNearby = 1;

        transformedRowInfo.forEach((d, i) => {
            // To aggregate bars, check if there is a same category on the next row.
            if(i + 1 < transformedRowInfo.length && d[field] === transformedRowInfo[i+1][field]) {
                if(aggregateStartIdx === -1) {
                    aggregateStartIdx = i;
                }
                sameCategoriesNearby++;
                return;
            }

            const barTop = aggregateStartIdx !== -1 ? yScale(aggregateStartIdx) : yScale(i);
            const barHeight = rowHeight * sameCategoriesNearby;
            const barWidth = barAreaWidth;
            const barLeft = (isLeft ? width - barWidth : 0);
            const textLeft = (isLeft ? width - barWidth - margin : barWidth + margin);
            const color = colorScale(d[field]);

            const rect = two.makeRect(barLeft, barTop, barWidth, barHeight);
            rect.fill = color;

            if(hoverValue && d[field] === hoverValue) {
                const hoverBgRectLeft = (isLeft ? barLeft - textAreaWidth : barLeft + barWidth)
                const hoverBgRect = two.makeRect(hoverBgRectLeft, barTop, textAreaWidth, barHeight);
                hoverBgRect.fill = "#EBEBEB";
            }

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

        drawVisTitle(titleText, { two, isLeft, width, height, titleSuffix });

        two.update();
        return two.teardown;
    });
    
    drawRegister("TrackRowInfoVisNominalBar", draw);

    // Context menu.
    function onContextMenu(e){
        e.preventDefault();
        
        const mouseViewportX = e.pageX;
        const mouseViewportY = e.pageY;

        PubSub.publish(EVENT.CONTEXT_MENU, {
            x: mouseViewportX,
            y: mouseViewportY,
            menuType: CONTEXT_MENU_TYPE.NOMINAL_BAR,
            field,
            value: hoverValue,
            onFilter: (f,t,v) => onFilterRows(f,t,v),
            onHighlight: (f,t,v) => onSearchRows(f,t,v)
        });    
    }

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
                fieldVal = transformedRowInfo[y][field];
                setHoverValue(fieldVal);
            } else {
                setMouseX(null);
                setHoverValue(null);
                destroyTooltip();
                return;
            }

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;
            
            PubSub.publish(EVENT.TOOLTIP, {
                x: mouseViewportX,
                y: mouseViewportY,
                content: <TooltipContent 
                    title={field}
                    value={fieldVal}
                    color={colorScale(fieldVal)}
                />
            });
        });

        // Handle mouse leave.
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => {
            setMouseX(null);
            setHoverValue(null);
        });

        // Clean up.
        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, transformedRowInfo, hoverValue]);

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
                onContextMenu={onContextMenu}
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
                onFilterRows={onFilterRows}
                transformedRowInfo={transformedRowInfo}
            />
        </div>
    );
}