import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./constants.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { drawVisTitle } from "./utils/vis.js";

import TrackControl from './TrackControl.js';

const margin = 5;

export default function TrackRowInfoVisLink(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        rowInfo,
        register,
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();
    const [mouseX, setMouseX] = useState(null);

    // Data, layouts and styles
    const { field, type, title } = fieldInfo;
    const isNominal = type === "nominal";
    
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
    
    register("TrackRowInfoVisBar", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const mouse = d3.mouse(canvas);
            const mouseX = mouse[0];
            const mouseY = mouse[1];

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
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => setMouseX(null));

        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height]);

    console.log("TrackRowInfoVis.render");
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
            <TrackControl
                isLeft={isLeft}
                isVisible={mouseX !== null}
                fieldInfo={fieldInfo}
            />
        </div>
    );
}