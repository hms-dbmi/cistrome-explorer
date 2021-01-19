import React, { useRef, useCallback, useEffect, useState, useMemo } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT, CONTEXT_MENU_TYPE } from "./utils/constants.js";
import { drawVisTitle } from "./utils/vis.js";

import TrackRowInfoControl from './TrackRowInfoControl.js';
import { TooltipContent, destroyTooltip } from "./Tooltip.js";
import { FILTER, HIGHLIGHTER, ARROW_UP, ARROW_DOWN } from './utils/icons.js';
import { getAggregatedValue } from "./utils/aggregate.js";
import { drawRowHighlightRect } from "./utils/linking.js";
import { HIGLASSMETA_DEFAULT } from "./utils/visualization-properties.js";

export const margin = 5;

/**
 * Component for visualization of row info quantitative or nominal attribute values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {number} titleHeight The height of the track title.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {boolean} isShowControlButtons Determine if control buttons should be shown.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {array} selectedRows The array of selected indices. 
 * @prop {array} highlitRows The array of highlit indices.
 * @prop {string} titleSuffix The suffix of a title, information about sorting and filtering status.
 * @prop {object} sortInfo The options for sorting rows of the field used in this track.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a highlight interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {boolean} helpActivated Whether to show help instructions or not.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisNominalBar(props) {
    const {
        left, top, width, height, titleHeight,
        field, type, alt, title, aggFunction, resolveYScale,
        isLeft,
        isShowControlButtons,
        rowInfo,
        transformedRowInfo,
        selectedRows,
        highlitRows,
        titleSuffix,
        sortInfo,
        filterInfo,
        onAddTrack,
        onSortRows,
        onHighlightRows,
        onFilterRows,
        helpActivated,
        drawRegister,
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();
    const [hoverValue, setHoverValue] = useState(null);

    // Data, layouts and styles
    const aggValue = d => getAggregatedValue(d, field, "nominal", aggFunction);

    const yScale = d3.scaleBand()
        .domain(range(transformedRowInfo.length))
        .range([titleHeight, height]);
    const rowHeight = yScale.bandwidth();

    const colorScale = useMemo(() => 
        d3.scaleOrdinal()
            .domain(Array.from(new Set(
                rowInfo.map(d => aggValue(d)))
            ).sort())
            .range([...d3.schemeTableau10]),
    [rowInfo]);
    // '#2299DB', '#FC5D5C',

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        drawRowHighlightRect(
            two, 
            selectedRows, 
            highlitRows, 
            titleHeight, 
            width, 
            height - titleHeight
        );
                
        const textAreaWidth = width - HIGLASSMETA_DEFAULT.TRACK.MIN_WIDTH;
        const showTextLabel = textAreaWidth > 0;
        const barAreaWidth = width - textAreaWidth;
        const fontSize = 12;

        // Render visual components for each row (i.e., bars and texts).
        const textAlign = isLeft ? "end" : "start";
        let aggregateStartIdx = -1, sameCategoriesNearby = 1;
        transformedRowInfo.forEach((d, i) => {
            const category = aggValue(d);

            // To aggregate bars, check if there is a same category on the next row.
            let isAggregateNext = false;
            if(
                i + 1 < transformedRowInfo.length
                && category === aggValue(transformedRowInfo[i+1])
            ) {
                if(aggregateStartIdx === -1) {
                    aggregateStartIdx = i;
                }
                sameCategoriesNearby++;
                isAggregateNext = true;
            }

            if(isAggregateNext) {
                return;
            }

            const barTop = (aggregateStartIdx !== -1 ? yScale(aggregateStartIdx) : yScale(i));
            const barHeight = rowHeight * sameCategoriesNearby;
            const barWidth = barAreaWidth;
            const barLeft = (isLeft ? width - barWidth : 0);
            const textLeft = (isLeft ? width - barWidth - margin : barWidth + margin);
            const color = colorScale(category);

            const rect = two.makeRect(barLeft, barTop, barWidth, barHeight);
            rect.fill = color;

            // Render text labels when the space is enough.
            if(barHeight >= fontSize && showTextLabel){
                const text = two.makeText(textLeft, barTop + barHeight/2, textAreaWidth, barHeight, category);
                text.fill = d3.hsl(color).darker(3);
                text.fontsize = fontSize;
                text.align = textAlign;
                text.baseline = "middle";
                text.overflow = "ellipsis";
            }

            aggregateStartIdx = -1;
            sameCategoriesNearby = 1;
        });

        // if(!isShowControlButtons) {
            drawVisTitle(title, { two, isLeft, width, height, titleSuffix });
        // }

        two.update();
        return two.teardown;
    });
    
    drawRegister("TrackRowInfoVisNominalBar", draw);

    // Context menu.
    function onContextMenu(e){
        e.preventDefault();
        
        const mouseViewportX = e.clientX;
        const mouseViewportY = e.clientY;

        const notOneOf = Array.from(colorScale.domain());
        notOneOf.splice(notOneOf.indexOf(hoverValue), 1);
        
        PubSub.publish(EVENT.CONTEXT_MENU, {
            x: mouseViewportX,
            y: mouseViewportY,
            menuType: CONTEXT_MENU_TYPE.NOMINAL_BAR,
            title: `Selected category: ${hoverValue}`,
            items: [
                // TODO: Temporally removing a highlighting feature
                // { title: "Highlight Rows", icon: HIGHLIGHTER, action: () => onHighlightRows(field, "nominal", hoverValue) },
                { title: "Filter Rows", icon: FILTER, action: () => onFilterRows(field, "nominal", notOneOf) },
                // { title: "Add Top Track with Selected Rows", icon: ARROW_UP, action: () => onAddTrack(field, "nominal", notOneOf, "top", hoverValue) },
                // { title: "Add Bottom Track with Selected Rows", icon: ARROW_DOWN, action: () => onAddTrack(field, "nominal", notOneOf, "bottom") }
            ]
        });    
    }

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            let hoveredCategory;
            if(y !== undefined){
                hoveredCategory = aggValue(transformedRowInfo[y]);
                setHoverValue(hoveredCategory);
                onHighlightRows(field, "nominal", hoveredCategory);
            } else {
                setHoverValue(null);
                destroyTooltip();
                onHighlightRows("");
                return;
            }

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;
            
            PubSub.publish(EVENT.TOOLTIP, {
                x: mouseViewportX,
                y: mouseViewportY,
                content: <TooltipContent 
                    title={title}
                    value={hoveredCategory}
                    color={colorScale(hoveredCategory)}
                />
            });
        });

        // Handle mouse leave.
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => {
            setHoverValue(null);
            onHighlightRows("");
        });

        // Clean up.
        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, transformedRowInfo, hoverValue, isShowControlButtons]);
    
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
                onContextMenu={onContextMenu}
                style={{
                    top: 0,
                    left: 0,
                    width: `${width}px`,
                    height: `${(height)}px`,
                    position: 'relative'
                }}
            />
            <TrackRowInfoControl
                isLeft={isLeft}
                isVisible={isShowControlButtons}
                field={field}
                type={type}
                title={title}
                aggFunction={aggFunction}
                top={titleHeight}
                searchLeft={left}
                sortAsceButtonHighlit={sortInfo && sortInfo.order === "ascending"}
                sortDescButtonHighlit={sortInfo && sortInfo.order === "descending"}
                filterButtonHighlit={filterInfo && filterInfo.notOneOf.length !== 0}
                onSortRows={onSortRows}
                onHighlightRows={onHighlightRows}
                onFilterRows={onFilterRows}
                filterInfo={filterInfo}
                transformedRowInfo={transformedRowInfo}
                rowInfo={rowInfo}
                helpActivated={helpActivated}
            />
        </div>
    );
}