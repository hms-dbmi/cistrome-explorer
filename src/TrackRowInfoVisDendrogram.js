import React, { useRef, useCallback, useEffect, useState } from "react";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { destroyTooltip } from "./utils/tooltip.js";
import { drawVisTitle } from "./utils/vis.js";
import { matrixToTree } from './utils/tree.js';
import { SORT_TREE } from './utils/icons.js';
import { TooltipContent } from './Tooltip.js';

import TrackRowInfoControl from './TrackRowInfoControl.js';

/**
 * Component for visualization of row info hierarchies.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} transformedRowInfo Array of JSON objects, one object for each row.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {function} onSortRows
 * @prop {function} onFilterRows
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisDendrogram(props) {
    const {
        left, top, width, height,
        fieldInfo,
        isLeft,
        transformedRowInfo,
        onSortRows,
        onFilterRows,
        drawRegister,
    } = props;
    
    // DOM refs.
    const divRef = useRef();
    const canvasRef = useRef();

    // Data refs.
    const descendantsRef = useRef();
    const delaunayRef = useRef();

    const [isMouseHover, setIsMouseHover] = useState(null);
    const [highlightNodeX, setHighlightNodeX] = useState(null);
    const [highlightNodeY, setHighlightNodeY] = useState(null);

    // Data, layouts and styles
    const { field } = fieldInfo;
    const isNominal = false;

    // Process the hierarchy data. Result will be null if the tree leaves
    // cannot be aligned based on the current rowInfo ordering.
    const hierarchyData = matrixToTree(transformedRowInfo.map(d => d[field]));
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

    function pointFromNode(d) {
        return [ (isLeft ? d.y : width - d.y), d.x ];
    }

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

        if(cannotAlign) {
            const rect = two.makeRect(0, 0, width, height);
            rect.fill = "white";
            rect.opacity = 0.8;
        }

        drawVisTitle(field, { two, isLeft, width, height });

        const points = descendants.map(pointFromNode);
        const delaunay = d3.delaunay.from(points);
        
        // Store delaunay object and descendants node array to references.
        descendantsRef.current = descendants;
        delaunayRef.current = delaunay;

        two.update();
        return two.teardown;
    });

    drawRegister("TrackRowInfoVisDendrogram", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            setIsMouseHover(true);

            if(cannotAlign) {
                // Show a tooltip indicating that the dendrogram has been hidden.
                const mouseViewportX = d3.event.clientX;
                const mouseViewportY = d3.event.clientY;
                PubSub.publish(EVENT.TOOLTIP, {
                    x: mouseViewportX,
                    y: mouseViewportY,
                    content: <TooltipContent 
                        title={`Dendrogram is inaccurate due to the current row ordering.`}
                    />
                });

                setHighlightNodeX(null);
                setHighlightNodeY(null);
            } else {
                // The dendrogram is visible, so no tooltip should be shown on hover.
                destroyTooltip();

                // Show hover indicator.
                const [mouseX, mouseY] = d3.mouse(canvas);
                const i = delaunayRef.current.find(mouseX, mouseY);
                const d = descendantsRef.current[i];
                const [pointX, pointY] = pointFromNode(d);

                setHighlightNodeX(pointX);
                setHighlightNodeY(pointY);
            }
        });

        d3.select(canvas).on("click", () => {
            if(!cannotAlign && delaunayRef.current && descendantsRef.current) {
                // Filter to select a subtree based on clicking near its root node.
                const [mouseX, mouseY] = d3.mouse(canvas);

                const i = delaunayRef.current.find(mouseX, mouseY);

                const subtree = [];
                let node = descendantsRef.current[i];
                while(node.parent) {
                    subtree.push(node.data.name);
                    node = node.parent;
                }
                subtree.reverse();
                onFilterRows(field, "tree", subtree);

                setHighlightNodeX(null);
                setHighlightNodeY(null);
            }
        });

        d3.select(canvas).on("mouseout", destroyTooltip);
        d3.select(div).on("mouseleave", () => setIsMouseHover(null));

        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [width, height, root]);

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
                    position: "absolute",
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`
                }}
            />
            {cannotAlign ? (
                <div onClick={() => onSortRows(fieldInfo.field, fieldInfo.type, "ascending")}
                    style={{
                        position: "absolute",
                        top: `${height / 2.0 - 17}px`,
                        left: `${width / 2.0 - 17}px`
                    }}>
                    <svg className={`chgw-button-alert`}
                        viewBox={SORT_TREE.viewBox}>
                        <title>{"Sort rows by hierarchy leaf order"}</title>
                        <path d={SORT_TREE.path} fill="currentColor"/>
                    </svg>
                </div>
            ) : null}
            {(isMouseHover && highlightNodeX && highlightNodeY) ? (
                <svg style={{
                    position: "absolute",
                    top: 0,
                    left: 0,
                    width: `${width}px`,
                    height: `${height}px`,
                    pointerEvents: "none"
                }}>
                    <circle 
                        r="6" 
                        fill="gray" 
                        stroke="none" 
                        opacity="0.6"
                        cx={highlightNodeX} 
                        cy={highlightNodeY}
                    />
                </svg>
            ) : null}
        </div>
    );
}
