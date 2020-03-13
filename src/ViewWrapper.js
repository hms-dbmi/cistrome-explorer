import React, { useRef, useState, useEffect, useCallback } from 'react';
import d3 from './utils/d3.js';
import uuidv4 from 'uuid/v4';
import ViewColumnBrush from './ViewColumnBrush.js';
import DataTableForIntervalTFs from './DataTableForIntervalTFs.js';
import './ViewWrapper.scss';
import { getRange } from './utils/viewport.js';

/**
 * Component for rendering genome interval selection tools.
 * @prop {object} viewBoundingBox The bounding box (i.e., left, top, width, height) of a HiGlass view.
 * @prop {object[]} viewportTracks An array of `viewport-projection-horizontal` track objects.
 * @prop {object} multivecTrack A object of `horizontal-multivec` track in the same view.
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval.
 * @prop {function} onViewportRemove The function to call upon removing a viewport track.
 * @prop {function} onRequestIntervalTFs The function to call upon making a request for further interval transcription factor data.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function ViewWrapper(props) {

    const {
        viewBoundingBox,
        viewportTracks,
        multivecTrack,
        onSelectGenomicInterval,
        onViewportRemove,
        onRequestIntervalTFs,
        drawRegister
    } = props;
    
    const divRef = useRef();
    const dragX = useRef(null);
    const newViewportId = useRef(null);
    const [mouseHoverX, setMouseHoverX] = useState(null);
    const [brushStartX, setBrushStartX] = useState(null);
    const [brushEndX, setBrushEndX] = useState(null);
    
    const [requestedIntervalParams, setRequestedIntervalParams] = useState(null);

    const brushBarHeight = 24;

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const div = divRef.current;
        const [mouseX, mouseY] = d3.mouse(div);

        dragX.current = mouseX;
        newViewportId.current = uuidv4();
        
        setBrushStartX(mouseX);
        setMouseHoverX(null);
    }, [dragX]);

    const dragged = useCallback(() => {
        const [mouseX, mouseY] = d3.mouse(divRef.current);

        const {start, end} = getRange(dragX.current, mouseX, 0, viewBoundingBox.width);
        
        setBrushStartX(start);
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

    // All hooks must be above this return statement, since they need to be executed in the same order.
    if(!viewBoundingBox || !multivecTrack || !multivecTrack.tilesetInfo) {
        // The view info has not yet loaded.
        return null;
    }

    const { top, left, width, height } = viewBoundingBox;
    const brushBarTop = height + 4;

    return (
        <div className="cistrome-hgw-view-wrapper">
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
                }}>
                {mouseHoverX ? 
                    <div className="col-tools-hover-line" style={{
                        left: `${mouseHoverX}px`
                    }}/>
                : null}
                {brushStartX !== null && brushEndX !== null ?
                    <div className="col-tools-brush" style={{
                        left: brushStartX,
                        width: brushEndX - brushStartX
                    }}/>
                : null}
                
                {viewportTracks ? 
                    (viewportTracks.map((viewportTrack, i) => {
                        if(!viewportTrack) {
                            return null;
                        }
                        return (
                            <ViewColumnBrush
                                key={i}
                                viewBoundingBox={viewBoundingBox}
                                viewportTrack={viewportTrack}
                                multivecTrack={multivecTrack}
                                onViewportRemove={onViewportRemove}
                                onRequestIntervalTFs={onRequestIntervalTFs}
                                onRequestIntervalTFs={(intervalParams) => {
                                    setRequestedIntervalParams(intervalParams);
                                }}
                            />
                        );
                    }))
                : null}
            </div>
            <div className="col-tools-hg-overlay">
                {mouseHoverX ? 
                    <div className="col-tools-hover-line" style={{
                        left: `${mouseHoverX}px`
                    }}/>
                : null}
                {brushStartX !== null && brushEndX !== null ?
                    <div className="col-tools-brush" style={{
                        left: brushStartX,
                        width: brushEndX - brushStartX
                    }}/>
                : null}
            </div>
            {requestedIntervalParams ? 
                <DataTableForIntervalTFs
                    left={left}
                    top={height + 40}
                    width={width}
                    height={600}
                    intervalParams={requestedIntervalParams}
                />
            : null}
        </div>
    );
};
