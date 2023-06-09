import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { drawVisTitle } from "./utils/vis.js";
import { TooltipContent, destroyTooltip } from "./Tooltip.jsx";
import TrackRowInfoControl from "./TrackRowInfoControl.jsx";
import { getAggregatedValue } from "./utils/aggregate.js";

const margin = 5;
// Color for the initial / undefined state.
const UNDEFINED_COLOR = "#dedede";

function createTooltipText(domain) {
    let text = "Click ";
    domain.forEach((domainVal, i) => {
        text += `${i+1}x to mark ${domainVal}`;
        if(i < domain.length - 1) {
            text += ", ";
        } else {
            text += ".";
        }
    });
    return text;
}

/**
 * Component for visualization of "dynamic" / user-defined nominal values, which can be changed by clicking.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {number} titleHeight The height of the track title.
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
        left, top, width, height, titleHeight,
        field, type, alt, title, aggFunction, resolveYScale,
        domain: nominalDomain, range: nominalRange,
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
    const [idToDomainValueMap, setIdToDomainValueMap] = useState({});

    // Range and domain length should match.
    console.assert(nominalRange.length === nominalDomain.length);

    function getRowColor(rowId) {
        const domainValue = idToDomainValueMap[rowId];
        if(domainValue !== undefined) {
            const domainIndex = nominalDomain.indexOf(domainValue);
            if(0 <= domainIndex && domainIndex < nominalRange.length) {
                const rangeValue = nominalRange[domainIndex];
                return rangeValue;
            }
        }
        return UNDEFINED_COLOR;
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
            const rowId = aggValue(info, field);

            const rowRect = two.makeRect(0, textTop, width, rowHeight);
            rowRect.fill = getRowColor(rowId);
            rowRect.stroke = null;
            rowRect.opacity = 0.7 + (hoverIndex !== null && i === hoverIndex ? 0.3 : 0);

            if(shouldRenderText && isTextLabel) {
                const text = two.makeText(textLeft, 30 + textTop + rowHeight/2, width, rowHeight, rowId);
                text.fill = "#333";
                text.fontsize = fontSize;
                text.align = textAlign;
                text.baseline = "middle";
                text.overflow = "ellipsis";
            }
        });
        // if(!isShowControlButtons) {
        drawVisTitle(title, { two, isLeft, width, height, titleSuffix });
        // }
        
        two.update();
        return two.teardown;
    });
    

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            let rowId;
            if(y !== undefined) {
                rowId = aggValue(transformedRowInfo[y], field);
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
                    title={createTooltipText(nominalDomain)}
                    value={rowId + (idToDomainValueMap[rowId] !== undefined ? ` (${idToDomainValueMap[rowId]})` : "")}
                />
            });
        });

        // Handle mouse click interaction to visit links.
        d3.select(canvas).on("click", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            if(y !== undefined) {
                const info = transformedRowInfo[y];
                const rowId = aggValue(info, field);
                setIdToDomainValueMap(prev => {
                    const next = Object.assign({}, prev);
                    if(next[rowId] === undefined) {
                        next[rowId] = (nominalDomain.length > 0 ? nominalDomain[0] : undefined);
                    } else {
                        const prevIndex = nominalDomain.indexOf(next[rowId]);
                        if(prevIndex === -1 || prevIndex >= nominalDomain.length - 1) {
                            next[rowId] = undefined;
                        } else {
                            next[rowId] = nominalDomain[prevIndex+1];
                        }
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
    }, [top, left, width, height, transformedRowInfo, hoverIndex, idToDomainValueMap]);

    return (
        <div
            ref={divRef}
            style={{
                top: `${top}px`,
                position: "relative",
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
                    position: "relative",
                    cursor: "pointer"
                }}
            />
            <TrackRowInfoControl
                isLeft={isLeft}
                top={titleHeight}
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