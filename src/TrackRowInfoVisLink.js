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
import { drawRowHighlightRect } from "./utils/linking.js";

const margin = 5;

/**
 * Component for visualization of row info URL values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {array} selectedRows The array of selected indices. 
 * @prop {array} highlitRows The array of highlit indices.
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
export default function TrackRowInfoVisLink(props) {
    const {
        left, top, width, height,
        fieldInfo,
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
        drawRegister,
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();
    const [hoverIndex, setHoverIndex] = useState(null);

    // Data, layouts and styles
    const { field, title, aggFunction } = fieldInfo;
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

        if(shouldRenderText) {
            // There is enough height to render the text elements.
            transformedRowInfo.forEach((info, i) => {
                const textTop = yScale(i);
                const textLeft = isLeft ? width - margin : margin;
                const titleField = title ? title : field;
                const diplayText = isTextLabel ? aggValue(info, titleField) : "Link";
                const text = two.makeText(textLeft, textTop + rowHeight/2, width, rowHeight, diplayText);
                text.fill = "#23527C";
                text.fontsize = fontSize;
                text.align = textAlign;
                text.baseline = "middle";
                text.overflow = "ellipsis";

                if(hoverIndex !== null && hoverIndex === i) {
                    // There is both a hovered link and enough height to render text, so also render an underline for the hovered link.
                    const textDimensions = two.measureText(text);
                    const textUnderlineLeft = isLeft ? textLeft - textDimensions.width : textLeft;
                    const textUnderlineTop = textTop + (rowHeight / 2) + (textDimensions.height / 2);
                    const textUnderline = two.makeLine(textUnderlineLeft, textUnderlineTop, textUnderlineLeft + textDimensions.width, textUnderlineTop);
                    textUnderline.stroke = "#23527C";
                    textUnderline.linewidth = 1;
                }
            });
        }
        
        drawRowHighlightRect(two, selectedRows, highlitRows, width, height);

        if(!isShowControlButtons) {
            drawVisTitle(field, { two, isLeft, width, height, titleSuffix });
        }
        
        two.update();
        return two.teardown;
    });
    
    drawRegister("TrackRowInfoVisLink", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            let fieldVal;
            if(y !== undefined) {
                fieldVal = aggValue(transformedRowInfo[y], field);
                setHoverIndex(y);
                onHighlightRows(field, "nominal", fieldVal);
            } else {
                setHoverIndex(null);
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
                    title={field}
                    value={fieldVal}
                />
            });
        });

        // Handle mouse click interaction to visit links.
        d3.select(canvas).on("click", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            if(y !== undefined) {
               window.open(transformedRowInfo[y][field]);
            }
        });

        // Handle mouse leave.
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => {
            setHoverIndex(null);
            onHighlightRows("");
        });

        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, transformedRowInfo, hoverIndex, isShowControlButtons]);

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
                fieldInfo={fieldInfo}
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