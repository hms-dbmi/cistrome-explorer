import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { drawVisTitle } from "./utils/vis.js";
import { TooltipContent, destroyTooltip } from "./Tooltip.js";
import TrackRowInfoControl from './TrackRowInfoControl.js';
import { getAggregatedValue } from "./utils/aggregate.js";

const margin = 5;
const numStates = 3; // undefined, negative, positive

/**
 * Component for visualization of row info URL values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {string} titleSuffix The suffix of a title, information about sorting and filtering status.
 * @prop {object} sortInfo The options for sorting rows of the field used in this track.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a highlight interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisComparison(props) {
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
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();
    const [hoverIndex, setHoverIndex] = useState(null);
    const [idToClicks, setIdToClicks] = useState({});

    function getRowColor(rowLabel) {
        if(idToClicks[rowLabel] === 2) {
            return "#00ff00";
        } else if(idToClicks[rowLabel] === 1) {
            return "#ff0000";
        }
        return "#dedede";
    }
    function getRowState(rowLabel) {
        if(idToClicks[rowLabel] === 2) {
            return "positive";
        } else if(idToClicks[rowLabel] === 1) {
            return "negative";
        }
        return null;
    }

    // Data, layouts and styles
    const minTrackWidth = 40;
    const isTextLabel = width > minTrackWidth;
    const aggValue = (d, f) => getAggregatedValue(d, f, "nominal", aggFunction);
    
    const fontSize = 10;
    const textAlign = isLeft ? "end" : "start";

    // Scales
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

        const shouldRenderText = (rowHeight >= fontSize);

        if(hoverIndex !== null) {
            // There is currently a hovered element, so render a background rect.
            const bgRect = two.makeRect(0, yScale(hoverIndex), width, rowHeight);
            bgRect.fill = "#EBEBEB";
        }

        transformedRowInfo.forEach((info, i) => {
            const textTop = yScale(i);
            const textLeft = isLeft ? width - margin : margin;
            const rowLabel = aggValue(info, field);

            const rowRect = two.makeRect(0, textTop, width, rowHeight);
            rowRect.fill = getRowColor(rowLabel);
            rowRect.stroke = null;
            rowRect.opacity = 0.7 + (hoverIndex !== null && i === hoverIndex ? 0.3 : 0);

            if(shouldRenderText && isTextLabel) {
                const text = two.makeText(textLeft, textTop + rowHeight/2, width, rowHeight, rowLabel);
                text.fill = "#333";
                text.fontsize = fontSize;
                text.align = textAlign;
                text.baseline = "middle";
                text.overflow = "ellipsis";
            }
        });
        
        drawVisTitle(title, { two, isLeft, width, height, titleSuffix });
        
        two.update();
        return two.teardown;
    });
    

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);
        console.log("drew")

        d3.select(canvas).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            let fieldVal;
            if(y !== undefined) {
                fieldVal = aggValue(transformedRowInfo[y], field);
                setHoverIndex(y);
            } else {
                setHoverIndex(null);
                destroyTooltip();
                return;
            }

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;
            
            PubSub.publish(EVENT.TOOLTIP, {
                x: mouseViewportX,
                y: mouseViewportY,
                content: <TooltipContent 
                    title={"Click once to mark negative, twice to mark positive."}
                    value={fieldVal + (idToClicks[fieldVal] ? ` (${getRowState(fieldVal)})` : '')}
                />
            });
        });

        // Handle mouse click interaction to visit links.
        d3.select(canvas).on("click", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            if(y !== undefined) {
               const info = transformedRowInfo[y];
               const rowLabel = aggValue(info, field);
                setIdToClicks(prev => {
                    const next = Object.assign({}, prev);
                    if(next[rowLabel] === undefined) {
                        next[rowLabel] = 1;
                    } else {
                        next[rowLabel] = (next[rowLabel] + 1) % numStates;
                    }
                    return next;
                });
               destroyTooltip();
            }
        });

        // Handle mouse leave.
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => setHoverIndex(null));

        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, transformedRowInfo, hoverIndex, idToClicks]);

    return (
        <div
            ref={divRef}
            style={{
                position: 'relative',
                width: `${width}px`,
                height: `${height}px`
            }}
        >
            <canvas
                ref={canvasRef}
                style={{
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    position: 'relative',
                    cursor: 'pointer'
                }}
            />
            <TrackRowInfoControl
                isLeft={isLeft}
                isVisible={isShowControlButtons}
                field={field}
                type={type}
                title={title}
                aggFunction={aggFunction}
                sortAsceButtonHighlit={sortInfo && sortInfo.order === "ascending"}
                sortDescButtonHighlit={sortInfo && sortInfo.order === "descending"}
                filterButtonHighlit={filterInfo && filterInfo.notOneOf.length !== 0}
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