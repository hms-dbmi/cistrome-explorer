import React, { useRef, useCallback, useEffect, useState, useMemo } from "react";
import PubSub from "pubsub-js";

import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { EVENT, CONTEXT_MENU_TYPE } from "./utils/constants.js";
import { drawVisTitle } from "./utils/vis.js";
import { matrixToTree } from './utils/tree.js';
import { SORT_TREE, HIGHLIGHTER } from './utils/icons.js';
import { TooltipContent, destroyTooltip } from './Tooltip.js';
import TrackRowInfoControl from "./TrackRowInfoControl.js";
import { FILTER } from './utils/icons.js';
import { drawRowHighlightRect } from "./utils/linking.js";

/**
 * Component for visualization of row info hierarchies.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {array} selectedRows The array of selected indices. 
 * @prop {array} highlitRows The array of highlit indices.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {boolean} isShowControlButtons Determine if control buttons should be shown.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a highlight interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisDendrogram(props) {
    const {
        left, top, width, height, titleHeight,
        field, type, alt, title, aggFunction, resolveYScale,
        rowInfo,
        transformedRowInfo,
        filterInfo,
        isLeft,
        isShowControlButtons,
        selectedRows,
        highlitRows,
        onAddTrack,
        onSortRows,
        onHighlightRows,
        onFilterRows,
        helpActivated,
        rowAggregated,
        drawRegister,
    } = props;
    
    // DOM refs.
    const divRef = useRef();
    const canvasRef = useRef();
    const axisRef = useRef();
    const minSimBarRef = useRef();

    // Data refs.
    const descendantsRef = useRef();
    const delaunayRef = useRef();
    const ancestor = useRef();
    const minSimilarity = useRef();

    const [maxDistance, setMaxDistance] = useState(transformedRowInfo?.[0]?.[field]?.[0]?.dist);
    const [highlightNodeX, setHighlightNodeX] = useState(null);
    const [highlightNodeY, setHighlightNodeY] = useState(null);
    const [minSimBarLeft, setMinSimBarLeft] = useState(initialMinSimBarLeft);
    const [showMinSimBar, setShowMinSimBar] = useState(false);

    // Data, layouts and styles
    const axisHeight = 30;
    const titleAreaWidth = 20;
    const visWidth = width - titleAreaWidth;

    const initialMinSimBarLeft = useMemo(() => {
        return isLeft ? 20 : width - 20
    }, [isLeft, width]);

    useEffect(() => {
        if(!filterInfo || filterInfo.length === 0) {
            // Upon reset button, reset the minSimilarity value to zero to reset fading out dendrograms.
            minSimilarity.current = undefined;
        }
    }, [filterInfo]);

    useEffect(() => {
        // Handling cases where the width changes by a `reziser`.
        let newLeft;
        if(isLeft) {
            newLeft = titleAreaWidth + (1 - minSimilarity.current / maxDistance) * visWidth;
        } else {
            newLeft = minSimilarity.current / maxDistance * visWidth;
        }
        setMinSimBarLeft(newLeft ? newLeft : initialMinSimBarLeft);
    }, [width, titleAreaWidth, visWidth, maxDistance, isLeft, minSimilarity]);

    // Process the hierarchy data. Result will be null if the tree leaves
    // cannot be aligned based on the current rowInfo ordering.
    const hierarchyData = matrixToTree(
        transformedRowInfo.map(d => d[field])
            // TODO: Remove below when we can aggregate `tree`
            .filter(d => d !== undefined)
    );
    const root = d3.hierarchy(hierarchyData);
    const leaves = root.leaves().map(l => l.data);

    useEffect(() => {
        let newMaxDist = filterInfo?.minSimilarity;
        if(!newMaxDist) {
            newMaxDist = transformedRowInfo?.[0]?.[field]?.[0]?.dist;
        }
        setMaxDistance(newMaxDist);
    }, [filterInfo, transformedRowInfo, field]);
    
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
        .size([height - titleHeight, visWidth])
        .separation(() => 1);
    treeLayout(root);

    function pointFromNode(d) {
        return [ (isLeft ? d.y + titleAreaWidth : visWidth - d.y), d.x ];
    }

    const dragX = useRef(null);

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const event = d3.event;
        dragX.current = event.sourceEvent.clientX;
    }, [dragX])

    const ended = useCallback(() => {
        dragX.current = null;
    }, [maxDistance, width, dragX, minSimBarLeft]);

    const dragged = useCallback(() => {
        const event = d3.event;
        const diff = event.sourceEvent.clientX - dragX.current;
        let newLeft = minSimBarLeft + diff;
        if(newLeft < 0) {
            newLeft = 0;
        } else if (newLeft > width) {
            newLeft = width;
        }
        minSimilarity.current = isLeft ? 
            maxDistance * (1 - (newLeft - titleAreaWidth) / visWidth) : 
            maxDistance * (newLeft / visWidth);
        setMinSimBarLeft(newLeft);
        
        // TODO: We want to uncomment the below line, but that is somehow removing `ancestors` in `filterInfo`.
        // onHighlightRows(field, "tree", minSimilarity.current);
    }, [maxDistance, width, dragX, minSimBarLeft]);

    // Detect drag events for repositioning minimum similarity bar.
    useEffect(() => {
        const minSimBar = minSimBarRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(minSimBar).call(drag);

        return () => d3.select(minSimBar).on(".drag", null);
    }, [minSimBarRef, started, dragged, ended]);

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });
        
        drawRowHighlightRect(two, selectedRows, highlitRows, titleHeight, width, height - titleHeight);

        // Draw the dendrogram.
        const descendants = root.descendants();

        let pathFunction;
        if(maxDistance) {
            // Since our dendrogram is rotated, *.y is indicating position along x-axis.
            if(isLeft){
                pathFunction = (d) => {
                    if(filterInfo?.minSimilarity && d.parent.data.dist > filterInfo.minSimilarity) return null;
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
                    if(filterInfo?.minSimilarity && d.parent.data.dist > filterInfo.minSimilarity) return null;
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
                        visWidth - d.parent.y, top + d.parent.x + titleHeight + 30,
                        visWidth - d.parent.y, top + d.x + titleHeight + 30,
                        visWidth - d.y, top + d.x + titleHeight + 30,
                        visWidth - d.parent.y, top + d.x + titleHeight + 30
                    );
                }
            }
        }
        descendants.forEach((d, i) => {
            if(i > 0) {
                const path = pathFunction(d);
                if(path) {
                    // A path that is out of a min similarity bar can be cut off.
                    path.stroke = "#555";
                    path.linewidth = 1;
                    path.opacity = showMinSimBar && d.parent?.data.dist > minSimilarity.current ? 0.2 : 1;
                }
            }
        });

        if(cannotAlign) {
            const rect = two.makeRect(0, 0, width, height);
            rect.fill = "#F6F6F6";
            rect.opacity = 1;
        }

        // if(!isShowControlButtons) {
            drawVisTitle(field, { two, isLeft, width, height });
        // }

        const points = descendants.map(pointFromNode);
        const delaunay = d3.delaunay.from(points);
        
        // Store delaunay object and descendants node array to references.
        descendantsRef.current = descendants;
        delaunayRef.current = delaunay;

        two.update();
        return two.teardown;
    });

    const drawAxis = useCallback((domElement) => {
        if(!maxDistance) return () => { };
        
        d3.select(domElement).selectAll("*").remove();

        const axisScale = d3.scaleLinear()
            .domain(isLeft ? [maxDistance, 0] : [0, maxDistance])
            .range([0, visWidth]);
        const axis = d3.axisBottom(axisScale)
            .ticks(Math.ceil(visWidth / 60))
            .tickFormat(d3.format("~s"));
        
        d3.select(domElement)
            .attr("width", visWidth)
            .attr("height", axisHeight + titleHeight)
            .append("g")
                .attr("transform", `translate(${-1}, ${titleHeight})`)
                .call(axis);
        
        d3.select(domElement)
            .selectAll("text")
                .attr("transform", `translate(${isLeft ? -3 : 3}, 0)`);

        return () => { /* Teardown */ };
    });

    drawRegister("TrackRowInfoVisDendrogram", draw);
    drawRegister("TrackRowInfoVisDendrogramAxis", drawAxis);

    // Context menu.
    function onContextMenuFromBranch(e){
        e.preventDefault();

        if(!cannotAlign && delaunayRef.current && descendantsRef.current && !showMinSimBar) {
            const mouseViewportX = e.clientX;
            const mouseViewportY = e.clientY;
            
            let node = ancestor.current;
            const ancestors = [];
            while(node.parent) {
                ancestors.push(node.data.name);
                node = node.parent;
            }
            ancestors.reverse();
            PubSub.publish(EVENT.CONTEXT_MENU, {
                x: mouseViewportX,
                y: mouseViewportY,
                title: "Options for dendrogram",
                menuType: CONTEXT_MENU_TYPE.TREE_ANCESTOR,
                items: [
                    { title: "Highlight Rows", icon: HIGHLIGHTER, action: () => onHighlightRows(field, "tree", ancestors) },
                    { title: "Filter Rows", icon: FILTER, action: () => onFilterRows(field, "tree", ancestors, false) }
                ]
            });
        }   
    }

    function onContextMenuFromMinSimBar(e){
        e.preventDefault();

        if(!cannotAlign && showMinSimBar) {
            const mouseViewportX = e.clientX;
            const mouseViewportY = e.clientY;
            
            PubSub.publish(EVENT.CONTEXT_MENU, {
                x: mouseViewportX,
                y: mouseViewportY,
                title: "Options for dendrogram",
                menuType: CONTEXT_MENU_TYPE.TREE_ANCESTOR,
                items: [
                    { title: "Highlight Rows", icon: HIGHLIGHTER, action: () => onHighlightRows(field, "tree", minSimilarity.current) },
                    { title: "Filter Rows", icon: FILTER, action: () => {
                        onFilterRows(field, "tree", minSimilarity.current, false);
                        setMinSimBarLeft(initialMinSimBarLeft); // Move to the original position.
                    }}
                ]
            });
        }   
    }

    useEffect(() => {
        const canvas = canvasRef.current;
        const svg = axisRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);
        const teardownSvg = drawAxis(svg);

        d3.select(canvas).on("mousemove", () => {

            const mouseViewportX = d3.event.clientX;
            const mouseViewportY = d3.event.clientY;

            if(cannotAlign) {
                setHighlightNodeX(null);
                setHighlightNodeY(null);
            } else {
                PubSub.publish(EVENT.TOOLTIP, {
                    x: mouseViewportX,
                    y: mouseViewportY,
                    content: <TooltipContent 
                        title={"Right-click to view highlight and filter options."}
                    />
                });

                // The dendrogram is visible, so show a tooltip with information about right-clicking.
                if(showMinSimBar) {
                    // Do not want to select nearest branch when min similarity bar shown.
                    setHighlightNodeX(null);
                    setHighlightNodeY(null);
                    return;
                }

                // Show hover indicator.
                const [mouseX, mouseY] = d3.mouse(canvas);
                const i = delaunayRef.current.find(mouseX, mouseY);
                const d = ancestor.current = descendantsRef.current[i];
                const [pointX, pointY] = pointFromNode(d);
                
                setHighlightNodeX(pointX);
                setHighlightNodeY(pointY);
                
                const ancestors = [];
                let node = d;
                while(node.parent) {
                    ancestors.push(node.data.name);
                    node = node.parent;
                }
                ancestors.reverse();
                onHighlightRows(field, "tree", ancestors);
            }
        });

        d3.select(canvas).on("mouseout", () => {
            destroyTooltip();
            if(!showMinSimBar) {
                onHighlightRows("");
            }
        });

        return () => {
            teardown();
            teardownSvg();
            d3.select(div).on("mouseleave", null);
        };
    }, [width, height, root, isShowControlButtons]);

    // Create the minimum similarity bar element.
    const minSimBar = useMemo(() => {
        return (
            <div
                ref={minSimBarRef}
                className="minimum-similarity-bar"
                onContextMenu={onContextMenuFromMinSimBar}
                style={{
                    visibility: showMinSimBar && !cannotAlign ? "visible" : "hidden",
                    top: "0px",
                    left: minSimBarLeft,
                    height: `${height}px`,
                    width: "2px"
                }}
            />
        );
    }, [minSimBarLeft, height, showMinSimBar, minSimBarRef, cannotAlign]);

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
                onContextMenu={onContextMenuFromBranch}
                style={{
                    position: "absolute",
                    top: 0,
                    width: `${width}px`,
                    height: `${height}px`
                }}
            />
            {cannotAlign ? null : 
                <TrackRowInfoControl
                    isLeft={isLeft}
                    top={titleHeight}
                    isVisible={isShowControlButtons}
                    field={field}
                    type={type}
                    title={title}
                    aggFunction={aggFunction}
                    searchLeft={left}
                    onFilterRows={onFilterRows}
                    rowInfo={rowInfo}
                    transformedRowInfo={transformedRowInfo}
                    filterButtonHighlit={showMinSimBar}
                    toggleMinSimBar={maxDistance && !cannotAlign ? () => {
                        // Show the minimum similarity bar only when similarity distance is available.
                        setMinSimBarLeft(initialMinSimBarLeft);
                        minSimilarity.current = undefined;
                        setShowMinSimBar(!showMinSimBar);
                    } : undefined}
                    helpActivated={helpActivated}
                />}
            {cannotAlign ? (
                <div style={{
                    position: "absolute",
                    top: `${30}px`,
                    left: '0px',
                    width: "100%",
                    border: '1px dotted gray',
                    height: "calc(100% - 30px)",
                    padding: "10px",
                    color: "gray",
                    fontWeight: "bold",
                    overflow: 'auto'
                }}>
                    Dendrogram is hidden since the ordering leaf nodes does not align with the current row ordering.
                    <br/>
                    <br/> 
                    💡
                    <span style={{fontStyle: 'italic', fontWeight: 'normal'}}>
                        {rowAggregated ? 
                            <>To see the dendrogram, you need to first deactivate <span style={{color: '#2299DB'}}>Aggregate By Tissue</span>.</> : 
                            <>You need to click on the button below to see the dendrogram.</>
                        }
                    </span>
                    
                    {!rowAggregated ?
                        <div onClick={() => onSortRows(field, type, "ascending")}
                            style={{
                                position: "absolute",
                                top: `${height / 2.0 - 17}px`,
                                left: `${width / 2.0 - 17}px`,
                                color: "lightgray"
                            }}>
                            <svg className={"hm-button-lg"}
                                viewBox={SORT_TREE.viewBox}>
                                <title>{"Sort rows by hierarchy leaf order"}</title>
                                <path d={SORT_TREE.path} fill="currentColor"/>
                            </svg>
                        </div> : null}
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
            {minSimBar}
        </div>
    );
}
