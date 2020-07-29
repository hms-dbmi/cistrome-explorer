import React, { useRef, useState, useEffect, useCallback } from 'react';
import d3 from './utils/d3.js';
import uuidv4 from 'uuid/v4';
import ViewColumnBrush from './ViewColumnBrush.js';
import { resolveIntervalCoordinates } from './utils/genome.js';
import { getRange } from './utils/viewport.js';
import { ARROW_H } from './utils/icons.js';
import './ViewWrapper.scss';

/**
 * Component for rendering genome interval selection tools.
 * @prop {object} viewBoundingBox The bounding box (i.e., left, top, width, height) of a HiGlass view.
 * @prop {object[]} viewportTracks An array of `viewport-projection-horizontal` track objects.
 * @prop {object} multivecTrack A object of `horizontal-multivec` track in the same view.
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval.
 * @prop {function} onViewportRemove The function to call upon removing a viewport track.
 * @prop {function} onGenomicIntervalSearch A function to call upon searching for TFs by using the selected interval. Optional.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function ViewWrapper(props) {

    const {
        viewBoundingBox,
        viewportTracks,
        multivecTrack,
        onSelectGenomicInterval,
        onViewportRemove,
        onGenomicIntervalSearch,
        onGenomicIntervalRowSort,
        drawRegister
    } = props;
    
    const divRef = useRef();
    const dragX = useRef(null);
    const newViewportId = useRef(null);
    const [assembly, setAssembly] = useState(null);
    const [mouseHoverX, setMouseHoverX] = useState(null);
    const [brushStartX, setBrushStartX] = useState(null);
    const [brushEndX, setBrushEndX] = useState(null);
    const [chrName, setChrName] = useState("");
    const [chrPos, setChrPos] = useState("");

    const brushBarHeight = 24;

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const div = divRef.current;
        const [mouseX, mouseY] = d3.mouse(div);

        dragX.current = mouseX;
        newViewportId.current = uuidv4();
        
        setBrushStartX(mouseX);
        setMouseHoverX(mouseX);
    }, [dragX]);

    const dragged = useCallback(() => {
        const [mouseX, mouseY] = d3.mouse(divRef.current);

        const {start, end} = getRange(dragX.current, mouseX, 0, viewBoundingBox.width);
        
        setBrushStartX(start);
        setMouseHoverX(mouseX);
        setBrushEndX(end);
    }, [dragX]);

    const ended = useCallback(() => {
        const [mouseX, mouseY] = d3.mouse(divRef.current);

        const { start, end } = getRange(dragX.current, mouseX, 0, viewBoundingBox.width);

        let startProp = start / viewBoundingBox.width;
        let endProp = end / viewBoundingBox.width;
        
        if(Math.abs(end - start) > 4) {
            onSelectGenomicInterval(startProp, endProp, newViewportId.current);
        }

        dragX.current = null;
        setBrushEndX(0);
        setBrushStartX(null);
    }, [dragX, brushStartX, brushEndX, viewBoundingBox]);

    useEffect(() => {
        const div = divRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);
        
        d3.select(div).call(drag);

        return () => d3.select(div).on(".drag", null);
    });

    useEffect(() => {
        const div = divRef.current;

        d3.select(div).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(div);
            setMouseHoverX(mouseX);
        });
        
        d3.select(div).on("mouseleave", () => {
            setMouseHoverX(null);
        });

        return () => {
            d3.select(div).on("mousemove", null);
            d3.select(div).on("mouseleave", null);
        };
    });

    useEffect(() => {
        try {
            setAssembly(multivecTrack.tilesetInfo.coordSystem);
        } catch(e) {
            console.log(e);   
        }
    });

    useEffect(() => {
        if(assembly && multivecTrack._xScale) {
            const domainX = multivecTrack._xScale.invert(mouseHoverX);
            resolveIntervalCoordinates(assembly, domainX)
                .then(result => {
                    setChrName(result[0][0]);
                    setChrPos(result[0][1].toLocaleString("en"));
                });
        }
    }, [mouseHoverX]);

    // All hooks must be above this return statement, since they need to be executed in the same order.
    if(!viewBoundingBox || !multivecTrack || !multivecTrack.tilesetInfo) {
        // The view info has not yet loaded.
        return null;
    }

    const { top, left, width, height } = viewBoundingBox;
    const brushBarTop = height + 4;

    return assembly ? (
        <div className="hm-view-wrapper">
            <div 
                style={{
                    position: "absolute",
                    top: `${top}px`,
                    left: `${left}px`,
                    width: `${width}px`, 
                    height: `${height}px`,
                    pointerEvents: "none"
                }}
            >
                <div className="col-tools-brush-bar" ref={divRef}
                    style={{
                        top: `${brushBarTop}px`,
                        height: `${brushBarHeight}px`,
                    }}
                >
                    {mouseHoverX && !brushStartX ? 
                        // Short vertical line.
                        <div className="col-tools-hover-line" 
                            style={{ left: `${mouseHoverX}px` }}
                        />
                    : null}
                    {brushStartX !== null && brushEndX !== null ?
                        // New narrow brushing rectangle.
                        <div className="col-tools-brush" 
                            style={{
                                left: brushStartX,
                                width: brushEndX - brushStartX
                            }}
                        />
                    : null}
                    {viewportTracks ? 
                        (viewportTracks.map((viewportTrack, i) => {
                            return viewportTrack ? (
                                // Narrow brushing rectangles based on viewports.
                                <ViewColumnBrush
                                    key={i}
                                    viewBoundingBox={viewBoundingBox}
                                    viewportTrack={viewportTrack}
                                    multivecTrack={multivecTrack}
                                    onViewportRemove={onViewportRemove}
                                    onGenomicIntervalSearch={onGenomicIntervalSearch}
                                    onGenomicIntervalRowSort={onGenomicIntervalRowSort}
                                />
                            ) : null;
                        }))
                    : null}
                    <span style={{
                        marginLeft: "4px",
                        verticalAlign: "middle",
                        color: "#b7b6b6",
                        display: "inline-block",
                        pointerEvents: "none"
                    }}>
                        <svg className={'hm-button-sm'}
                            style={{ color: "#b7b6b6", verticalAlign: "middle" }}
                            viewBox={ARROW_H.viewBox}>
                            <path d={ARROW_H.path} fill="currentColor"/>
                        </svg>
                        Genomic Interval Selection
                    </span> 
                </div>
                <div className="col-tools-hg-overlay">
                    {brushStartX !== null && brushEndX !== null ?
                        // New large brushing rectangle.
                        <div className="col-tools-brush-black" 
                            style={{
                                left: brushStartX,
                                width: brushEndX - brushStartX
                            }}
                        />
                    : null}
                    {mouseHoverX && brushStartX == null ? 
                        // Long vertical line.
                        <div className="col-tools-hover-line" 
                            style={{ left: `${mouseHoverX}px` }}
                        />
                    : null}
                    {mouseHoverX  ? 
                        // Text label of chromosome position.
                        <div className="col-tools-hover-line-info" 
                            style={{ top: "1px", left: `${mouseHoverX}px` }}
                        >
                            {`${chrName}: ${chrPos}`}
                        </div>
                    : null}
                </div>
            </div>
        </div>
    ) : null;
};
