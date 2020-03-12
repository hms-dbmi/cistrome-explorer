// TODO: Remove unused.
import React, { useRef, useState, useEffect, useCallback } from 'react';
import d3 from './utils/d3.js';
import uuidv4 from 'uuid/v4';

import './ViewColumnBrushTools.scss';

/**
 * Component for rendering genome interval selection tools.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {(object|null)} combinedTrack The `combined` track, parent of the `horizontal-multivec` track.
 * @prop {object[]} siblingTracks An array of `viewport-projection-horizontal` track objects, which
 *                                are siblings of `multivecTrack` (children of the same `combined` track).
 * @prop {string} colToolsPosition The value of the `colToolsPosition` option.
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval.
 * @prop {function} onRequestIntervalTFs The function to call upon making a request for further interval transcription factor data.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
// TODO: Update parameters 
export default function ViewColumnBrushTools(props) {

    const {
        viewBoundingBox,
        // trackAssembly,
        combinedTrack,
        siblingTracks,
        colToolsPosition,
        onSelectGenomicInterval,
        onRequestIntervalTFs,
        drawRegister
    } = props;

    const divRef = useRef();
    const dragX = useRef(null);
    const brushUid = useRef(null);
    const [mouseHoverX, setMouseHoverX] = useState(null);
    const [brushStartX, setBrushStartX] = useState(null);
    const [brushWidth, setBrushWidth] = useState(0);
    const [mouseChrX, setMouseChrX] = useState(null);
    
    // if(!trackAssembly) {
    //     console.warn("trackAssembly is null");
    //     return null;
    // }

    // TODO: Remove this.
    const isTop = (colToolsPosition === "top");

    const brushBarHeight = 24;

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const div = divRef.current;
        const [mouseX, mouseY] = d3.mouse(div);

        dragX.current = mouseX;
        brushUid.current = uuidv4();
        
        setBrushStartX(mouseX);
        setMouseHoverX(null);
    }, [dragX]);

    const dragged = useCallback(() => {
        const [mouseX, mouseY] = d3.mouse(divRef.current);
        const bWidth = mouseX - dragX.current;

        setBrushWidth(bWidth);
    }, [dragX]);

    const ended = useCallback(() => {
        const [mouseX, mouseY] = d3.mouse(divRef.current);

        const startProp = dragX.current / viewBoundingBox.width;
        const endProp = mouseX / viewBoundingBox.width;
        onSelectGenomicInterval(startProp, endProp, brushUid.current);

        dragX.current = null;
        setBrushWidth(0);
        setBrushStartX(null);
    }, [dragX, brushStartX, brushWidth, viewBoundingBox]);

    useEffect(() => {
        const div = divRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);
        
        d3.select(div).call(drag);

        return () => d3.select(div).on(".drag", null);
    }, [dragX, divRef]);

    useEffect(() => {
        const div = divRef.current;

        d3.select(div).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(div);

            // TODO: Change this to the actual chr position.
            // const xAbsPos = multivecTrack._xScale.invert(mouseX);

            setMouseHoverX(mouseX);
            // setMouseChrX(xAbsPos);
        });
        
        d3.select(div).on("mouseleave", () => {
            setMouseHoverX(null);
        });

        return () => {
            d3.select(div).on("mousemove", null);
            d3.select(div).on("mouseleave", null);
        };
    }, [dragX, divRef, mouseHoverX]);

    // All hooks must be above this return statement, since they need to be executed in the same order.
    if(!viewBoundingBox) {
        // The view info has not yet loaded.
        return null;
    }

    const { top, left, width, height } = viewBoundingBox;
    const brushBarTop = height + 4;

    return (
        <div style={{
            position: "absolute",
            top: `${top}px`,
            left: `${left}px`,
            width: `${width}px`, 
            height: `${height}px`,
            pointerEvents: "none"
        }}>
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
                {brushStartX && brushWidth !== 0 ?
                    <div className="col-tools-brush" style={{
                        left: brushStartX,
                        width: brushWidth
                    }}/>
                    : null}
                {mouseChrX ? 
                    <div className="col-tools-chr-info" style={{
                        left: `${mouseHoverX + 16}px`,
                        lineHeight: `${brushBarHeight}px`
                    }}>
                        {(+mouseChrX.toFixed()).toLocaleString("en")}
                    </div>
                    : null}
                {/* TODO: Remove Code Below when ready to release */}
                {/* <button 
                    className="col-tools-target"
                    onClick={() => onSelectGenomicInterval(0.25, 0.75)}
                    style={{
                        marginTop: (isTop ? (2 * brushBarHeight / 3) : 2)
                    }}
                >Select current interval</button>
                {siblingTracks ? (
                    <div className="col-tools-selection-info">
                        {siblingTracks.map((siblingTrack, i) => (
                            <TrackColSelectionInfo
                                key={i}
                                width={width}
                                height={brushBarHeight}
                                colToolsPosition={colToolsPosition}
                                // trackAssembly={trackAssembly}
                                projectionTrack={siblingTrack}
                                onRequestIntervalTFs={onRequestIntervalTFs}
                                drawRegister={drawRegister}
                            />
                        ))}
                    </div>
                ) : null} */}
                {/* TODO: remove above. */}
            </div>
            <div className="col-tools-hg-overlay">
                {mouseHoverX ? 
                    <div className="col-tools-hover-line" style={{
                        left: `${mouseHoverX}px`
                    }}/>
                    : null}
                {brushStartX && brushWidth !== 0 ?
                    <div className="col-tools-brush" style={{
                        left: brushStartX,
                        width: brushWidth
                    }}/>
                    : null}
            </div>
        </div>
    );
};