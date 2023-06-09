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
import { drawRowHighlightRect } from "./utils/linking.js";

const margin = 5;

/**
 * Component for visualization of row info URL values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {number} titleHeight The height of the track title.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {array} selectedRows The array of selected indices. 
 * @prop {array} highlitRows The array of highlit indices.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {string} shortName Alternative string value to show when the track is too narrow.
 * @prop {function} addTrackOnClick A callback function upon click on a row that adds a separate track on the top.
 * @prop {boolean} isLeft Is this view on the left side of the track?
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
export default function TrackRowInfoVisLink(props) {
    const {
        left, top, width, height, titleHeight,
        field, type, alt, title, aggFunction, resolveYScale, shortName = "Link", addTrackOnClick,
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
    const [hoverIndex, setHoverIndex] = useState(null);

    // Data, layouts and styles
    const minTrackWidth = 40;
    const isTextLabel = width > minTrackWidth;
    const aggValue = (d, f) => getAggregatedValue(d, f, "nominal", aggFunction);
    
    const fontSize = 10;
    const textAlign = isLeft ? "end" : "start";

    // Scales
    const yScale = d3.scaleBand()
        .domain(range(transformedRowInfo.length))
        .range([titleHeight, height]);
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
                const diplayText = isTextLabel ? aggValue(info, alt ? alt : field) : shortName;
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
        
        drawRowHighlightRect(
            two, 
            selectedRows, 
            highlitRows, 
            titleHeight, 
            width, 
            height - titleHeight
        );

        // if(!isShowControlButtons) {
        drawVisTitle(title, { two, isLeft, width, height, titleSuffix });
        // }
        
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
                    title={title}
                    value={fieldVal}
                />
            });
        });

        // Handle mouse click interaction to visit links.
        d3.select(canvas).on("click", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            const hoverValue = aggValue(transformedRowInfo[y], field);

            if(y !== undefined) {
                if(addTrackOnClick) {
                    const notOneOf = transformedRowInfo.map(d => aggValue(d, field));
                    notOneOf.splice(notOneOf.indexOf(hoverValue), 1);
                    onAddTrack(field, "nominal", notOneOf, "top", hoverValue);
                    console.log();
                }
                else {
                    window.open(hoverValue);
                }
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
                helpActivated={helpActivated}
            />
        </div>
    );
}