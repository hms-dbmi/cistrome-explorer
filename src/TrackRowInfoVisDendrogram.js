import React, { useRef, useCallback, useEffect, useState, useMemo } from "react";
import range from "lodash/range";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { drawVisTitle } from "./utils/vis.js";
import { matrixToTree } from './utils/tree.js';
import { EXCLAMATION } from './utils/icons.js';

import TrackRowInfoControl from './TrackRowInfoControl.js';

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
        onSortRows,
        drawRegister,
    } = props;
    
    const divRef = useRef();
    const canvasRef = useRef();

    const [mouseX, setMouseX] = useState(null);

    // Data, layouts and styles
    const { field } = fieldInfo;
    const isNominal = false;

    const yScale = d3.scaleBand()
        .domain(range(rowInfo.length))
        .range([0, height]);

    // Process the hierarchy data. Result will be null if the tree leaves
    // cannot be aligned based on the current rowInfo ordering.
    // const root = useMemo(() => {
    const hierarchyData = matrixToTree(rowInfo.map(d => d[field]));
    const root = d3.hierarchy(hierarchyData);
    const leaves = root.leaves().map(l => l.data);

    // Check whether dendrogram leaves can be aligned to the current row ordering.
    let cannotAlign = false;
    for(let i = 0; i < leaves.length; i++) {
        if(leaves[i].i !== i) {
            cannotAlign = true;
            break;
        }
    }
    const treeLayout = d3.cluster()
        .size([height, width])
        .separation(() => 1);
    treeLayout(root);

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });
        
        // Draw the dendrogram.
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

        if(cannotAlign){
            const rect = two.makeRect(0, 0, width, height);
            rect.fill = "white";
            rect.opacity = 0.8;
        }

        drawVisTitle(field, { two, isLeft, isNominal, width });

        two.update();
        return two.teardown;
    });

    drawRegister("TrackRowInfoVisDendrogram", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            setMouseX(true);

            if(cannotAlign) {
                // Show a tooltip indicating that the dendrogram has been hidden.
                const mouseViewportX = d3.event.clientX;
                const mouseViewportY = d3.event.clientY;
                PubSub.publish(EVENT.TOOLTIP, {
                    x: mouseViewportX,
                    y: mouseViewportY,
                    content: `Dendrogram is hidden due to the current row ordering.`
                });
            } else {
                // The dendrogram is visible, so no tooltip should be shown on hover.
                destroyTooltip();
            }
        });
        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => setMouseX(null));

        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [top, left, width, height, root]);

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
                searchTop={null}
                searchLeft={null}
                onSortRows={onSortRows}
                onSearchRows={null}
            />
            {cannotAlign ? 
                <div
                    style={{
                        position: "absolute",
                        pointerEvents: "none",
                        top: `${height / 2.0 - 25}px`,
                        left: `${width / 2.0 - 25}px`,
                    }}
                >
                    <svg 
                        className={`chgw-button-alert`}
                        viewBox={EXCLAMATION.viewBox}
                    >
                        <path d={EXCLAMATION.path} fill="currentColor"/>
                    </svg>
                </div>
                : null}
        </div>
    );
}