import React, { useRef, useCallback, useEffect, useState } from "react";
import Two from "./utils/two.js";
import range from "lodash/range";
import PubSub from "pubsub-js";
import d3 from "./utils/d3.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { EVENT } from "./constants.js";

import TrackControl from './TrackControl.js';
import TrackRowInfoVisBar from './TrackRowInfoVisBar.js';
import TrackRowInfoVisLink from './TrackRowInfoVisLink.js';
import TrackRowInfoVisDendrogram from './TrackRowInfoVisDendrogram.js';

const visTypeToVisComponent = {
    "nominal": TrackRowInfoVisBar,
    "quantitative": TrackRowInfoVisBar,
    "url": TrackRowInfoVisLink,
    "tree": TrackRowInfoVisDendrogram
};

export const margin = 5;

/**
 * Visualization for rendering an attribute based on its type.
 * @prop {object} two An instance of the Two class.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 */
export default function TrackRowInfoVis(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        rowInfo,
        register,
    } = props;

    const canvasRef = useRef();
    const [mouseX, setMouseX] = useState(null);

    // Data, layouts and styles
    const { field, type } = fieldInfo;
    const isNominal = type === "nominal";
    const barAreaWidth = isNominal ? 20 : width - 20;
    const titleFontSize = 12;

    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        // Title
        const titleLeft = (isLeft ? margin : width - margin);
        const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;
        const titleText = field;
        const title = two.makeText(titleLeft, top, 12, barAreaWidth, titleText);
        title.fill = "#9A9A9A";
        title.fontsize = titleFontSize;
        title.align = isLeft ? "end" : "start";
        title.baseline = "top";
        title.rotation = titleRotate;

        two.update();
        return two.teardown;
    });
    
    register("TrackRowInfoVis", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
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

        return teardown;
    }, [top, left, width, height]);

    console.log("TrackRowInfoVis.render");
    return (
        <div
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