import React, { useRef, useCallback, useEffect, useState } from "react";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT } from "./utils/constants.js";
import { drawVisTitle } from "./utils/vis.js";
import { matrixToTree } from './utils/tree.js';
import { SORT_TREE } from './utils/icons.js';
import { TooltipContent, destroyTooltip } from './Tooltip.js';
import TrackRowInfoControl from "./TrackRowInfoControl.js";

/**
 * Component for visualization of row info hierarchies.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after transforming by filtering and sorting according to the selected rows.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {boolean} isShowControlButtons Determine if control buttons should be shown.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisDendrogram(props) {
    const {
        left, top, width, height,
        transformedRowInfo,
        fieldInfo,
        isLeft,
        isShowControlButtons,
        onAddTrack,
        onSortRows,
        onFilterRows,
        drawRegister,
    } = props;
    
    // DOM refs.
    const divRef = useRef();
    const canvasRef = useRef();
    const axisRef = useRef();

    // Data refs.
    const descendantsRef = useRef();
    const delaunayRef = useRef();

    const [highlightNodeX, setHighlightNodeX] = useState(null);
    const [highlightNodeY, setHighlightNodeY] = useState(null);

    // Data, layouts and styles
    const { field } = fieldInfo;
    const axisHeight = 30;
    const titleAreaWidth = 20;
    const visWidth = width - titleAreaWidth;

    // Process the hierarchy data. Result will be null if the tree leaves
    // cannot be aligned based on the current rowInfo ordering.
    const hierarchyData = matrixToTree(transformedRowInfo.map(d => d[field]));
    const root = d3.hierarchy(hierarchyData);
    const leaves = root.leaves().map(l => l.data);
    const encodeDistance = hierarchyData.children[0] && hierarchyData.children[0].dist !== undefined;
    const maxDistance = encodeDistance ? hierarchyData.children[0].dist : null;

    // When subtrees are filtered, make sure the leaves touch the baseline.
    leaves.forEach(d => d.dist = 0);
    
    // Check whether dendrogram leaves can be aligned to the current row ordering.
    let cannotAlign = false;
    for(let i = 0; i < leaves.length; i++) {
        if(leaves[i].i !== i) {
            cannotAlign = true;
            break;
        }
    }
    const treeLayout = d3.cluster()
        .size([height, visWidth])
        .separation(() => 1);
    treeLayout(root);

    function pointFromNode(d) {
        return [ (isLeft ? d.y + titleAreaWidth : visWidth - d.y), d.x ];
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
        if(encodeDistance) {
            // Since our dendrogram is rotated, *.y is indicating position along x-axis.
            if(isLeft){
                pathFunction = (d) => {
                    d.y = visWidth * (1 - d.data.dist / maxDistance);
                    d.parent.y = visWidth * (1 - d.parent.data.dist / maxDistance);
                    return two.makePath(
                        titleAreaWidth + d.parent.y, top + d.parent.x,
                        titleAreaWidth + d.parent.y, top + d.x,
                        titleAreaWidth + d.y, top + d.x,
                        titleAreaWidth + d.parent.y, top + d.x
                    );
                }
            } else {
                pathFunction = (d) => {
                    d.y = visWidth * (1 - d.data.dist / maxDistance);
                    d.parent.y = visWidth * (1 - d.parent.data.dist / maxDistance);
                    return two.makePath(
                        visWidth - d.parent.y, top + d.parent.x,
                        visWidth - d.parent.y, top + d.x,
                        visWidth - d.y, top + d.x,
                        visWidth - d.parent.y, top + d.x
                    );
                }
            }
        } else {
            // TODO: Remove this part when we always encode similarity distance in dendrogram.
            if(isLeft){
                pathFunction = (d) => {
                    return two.makePath(
                        titleAreaWidth + d.parent.y, top + d.parent.x,
                        titleAreaWidth + d.parent.y, top + d.x,
                        titleAreaWidth + d.y, top + d.x,
                        titleAreaWidth + d.parent.y, top + d.x
                    );
                }
            } else {
                pathFunction = (d) => {
                    return two.makePath(
                        visWidth - d.parent.y, top + d.parent.x,
                        visWidth - d.parent.y, top + d.x,
                        visWidth - d.y, top + d.x,
                        visWidth - d.parent.y, top + d.x
                    );
                }
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

    const drawAxis = useCallback((domElement) => {
        if(!encodeDistance) return () => { };
        
        d3.select(domElement).selectAll("*").remove();

        const axisScale = d3.scaleLinear()
            .domain(isLeft ? [maxDistance, 0] : [0, maxDistance])
            .range([0, visWidth]);
        const axis = d3.axisBottom(axisScale)
            .ticks(Math.ceil(visWidth / 60))
            .tickFormat(d3.format("~s"));
        
        d3.select(domElement)
            .attr("width", visWidth)
            .attr("height", axisHeight)
            .append("g")
                .attr("transform", `translate(${-1}, 0)`)
                .call(axis);
        
        d3.select(domElement)
            .selectAll("text")
                .attr("transform", `translate(${isLeft ? -3 : 3}, 0)`);

        return () => { /* Teardown */ };
    });

    drawRegister("TrackRowInfoVisDendrogram", draw);
    drawRegister("TrackRowInfoVisDendrogramAxis", drawAxis);

    useEffect(() => {
        const canvas = canvasRef.current;
        const svg = axisRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);
        const teardownSvg = drawAxis(svg);

        d3.select(canvas).on("mousemove", () => {

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
                onFilterRows(field, "tree", subtree, false);

                setHighlightNodeX(null);
                setHighlightNodeY(null);
            }
        });

        d3.select(canvas).on("mouseout", destroyTooltip);

        return () => {
            teardown();
            teardownSvg();
            d3.select(div).on("mouseleave", null);
        };
    }, [width, height, root]);

    return (
        <div
            ref={divRef}
            style={{
                top: `${top}px`,
                position: 'relative',
                width: `${width}px`,
                height: `${height}px`,
            }}
        >
            <svg ref={axisRef} 
                style={{
                    left: isLeft ? titleAreaWidth : 0,
                    pointerEvents: "none",
                    position: "absolute"
                }}
            />
            <canvas
                ref={canvasRef}
                style={{
                    position: "absolute",
                    top: 0,
                    width: `${width}px`,
                    height: `${height}px`
                }}
            />
            <TrackRowInfoControl
                isLeft={isLeft}
                isVisible={isShowControlButtons}
                fieldInfo={fieldInfo}
                searchTop={top}
                searchLeft={left}
                onFilterRows={onFilterRows}
                transformedRowInfo={transformedRowInfo}
            />
            {cannotAlign ? (
                <div onClick={() => onSortRows(fieldInfo.field, fieldInfo.type, "ascending")}
                    style={{
                        position: "absolute",
                        top: `${height / 2.0 - 17}px`,
                        left: `${width / 2.0 - 17}px`,
                        color: "lightgray"
                    }}>
                    <svg className={"chw-button-lg"}
                        viewBox={SORT_TREE.viewBox}>
                        <title>{"Sort rows by hierarchy leaf order"}</title>
                        <path d={SORT_TREE.path} fill="currentColor"/>
                    </svg>
                </div>
            ) : null}
            {(isShowControlButtons && highlightNodeX && highlightNodeY) ? (
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
