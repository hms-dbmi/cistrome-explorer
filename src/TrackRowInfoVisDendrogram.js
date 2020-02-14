import React, { useRef, useCallback, useEffect, useState } from "react";
import range from "lodash/range";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { drawVisTitle } from "./utils/vis.js";
import { matrixToTree } from './utils/tree.js';

/**
 * Component for visualization of row info hierarchies.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {string} viewId The viewId for the horizontal-multivec track.
 * @prop {string} trackId The trackId for the horizontal-multivec track.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisDendrogram(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        viewId,
        trackId,
        rowInfo,
        drawRegister,
    } = props;

    // TODO: fix infinite recursion bug
   

    const divRef = useRef();
    const canvasRef = useRef();
    const [mouseX, setMouseX] = useState(null);

    // Data, layouts and styles
    const { field } = fieldInfo;
    const isNominal = false;

    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        drawVisTitle(field, { two, isLeft, isNominal, width });

        const hierarchyData = matrixToTree(rowInfo.map(d => d[field]));

        const root = d3.hierarchy(hierarchyData);

        const treeLayout = d3.cluster()
            .size([height, width])
            .separation(() => 1);
        treeLayout(root);

        const descendants = root.descendants();

        let pathFunction;
        if(isLeft){
            pathFunction = (d) => {
                return two.makePath(
                    d.parent.y, top + d.parent.x,
                    d.parent.y, top + d.x,
                    d.y, top + d.x,
                    d.parent.y, top + d.x
                );
            }
        } else {
            pathFunction = (d) => {
                return two.makePath(
                    width -  d.parent.y, top + d.parent.x,
                    width - d.parent.y, top + d.x,
                    width - d.y, top + d.x,
                    width - d.parent.y, top + d.x
                );
            }
        }
        
        descendants.forEach((d, i) => {
            if(i > 0) {
                const path = pathFunction(d);
                path.stroke = "#555";
                path.opacity = 0.6;
                path.linewidth = 1.5;
            }
        });
        
        two.update();
        return two.teardown;
    });
    
    drawRegister("TrackRowInfoVisDendrogram", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(canvas);
            const y = yScale.invert(mouseY);
            if(y !== undefined){
                setMouseX(true);
            } else {
                setMouseX(null);
                return;
            }
        });
        d3.select(div).on("mouseleave", () => setMouseX(null));

        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, rowInfo]);

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
        </div>
    );
}