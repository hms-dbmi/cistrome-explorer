import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./constants.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { drawVisTitle } from "./utils/vis.js";

import TrackRowInfoControl from './TrackRowInfoControl.js';

const margin = 5;

/**
 * Component for visualization of row info URL values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisLink(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        rowInfo,
        drawRegister,
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();
    const [mouseX, setMouseX] = useState(null);

    // Data, layouts and styles
    const { field, title } = fieldInfo;
    const isNominal = false;
    
    const fontSize = 10;
    const textAlign = isLeft ? "end" : "start";

    // Scales
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

        drawVisTitle(field, { two, isLeft, isNominal, width });

        if(rowHeight < fontSize) {
            // Bail out if there is not enough height per row to render links.
            return two.teardown;
        }

        rowInfo.forEach((info, i) => {
            const textTop = yScale(i);
            const textLeft = isLeft ? width - margin : margin;
            const titleField = title ? title : field;
    
            const text = two.makeText(textLeft, textTop + rowHeight/2, width, rowHeight, info[titleField]);
            text.fill = "#23527C";
            text.fontsize = fontSize;
            text.align = textAlign;
            text.baseline = "middle";
            text.overflow = "ellipsis";
        });
        
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

        // Handle mouse click interaction to visit links.
        d3.select(canvas).on("click", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);

            const y = yScale.invert(mouseY);
            if(y !== undefined) {
               window.open(rowInfo[y][field]);
            }
        });

        // Handle mouse leave.
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => setMouseX(null));

        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, rowInfo]);

    console.log("TrackRowInfoVisLink.render");
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
            />
        </div>
    );
}