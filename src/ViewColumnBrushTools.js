// TODO: Remove unused.
import React, { useRef, useState, useEffect, useCallback } from 'react';
import d3 from './utils/d3.js';
import uuidv4 from 'uuid/v4';
import { CLOSE, SEARCH } from './utils/icons.js';
import './ViewColumnBrushTools.scss';

/**
 * Component for rendering genome interval selection tools.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {(object|null)} combinedTrack The `combined` track, parent of the `horizontal-multivec` track.
 * @prop {object[]} viewportTracks An array of `viewport-projection-horizontal` track objects, which
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
        viewportTracks,
        colToolsPosition,
        onSelectGenomicInterval,
        onViewportRemove,
        onRequestIntervalTFs,
        drawRegister
    } = props;
    
    const divRef = useRef();
    const dragX = useRef(null);
    const brushUid = useRef(null);
    const [mouseHoverX, setMouseHoverX] = useState(null);
    const [brushStartX, setBrushStartX] = useState(null);
    const [brushWidth, setBrushWidth] = useState(0);

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

        // Handle the oposite direction of dragging along x-axis.
        setBrushStartX(bWidth < 0 ? mouseX : dragX.current);
        setBrushWidth(bWidth < 0 ? -bWidth : bWidth);
    }, [dragX]);

    const ended = useCallback(() => {
        const [mouseX, mouseY] = d3.mouse(divRef.current);

        let startProp = dragX.current / viewBoundingBox.width;
        let endProp = mouseX / viewBoundingBox.width;
        if(startProp > endProp) {
            const temp = startProp;
            startProp = endProp;
            endProp = temp;
        }
        if(startProp !== endProp) {
            onSelectGenomicInterval(startProp, endProp, brushUid.current);
        }

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
            setMouseHoverX(mouseX);
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
                
                {/* brushes for viewport-projection-horizontal tracks */}
                {viewportTracks ? 
                    (viewportTracks.map((viewportTrack, i) => {
                        if(!viewportTrack) return null;

                        const absDomain = viewportTrack.viewportXDomain;
                        const startX = viewportTrack._xScale(absDomain[0]);
                        const endX = viewportTrack._xScale(absDomain[1]);

                        let assembly = null;
                        try{
                            console.log(viewportTrack);
                            assembly = viewportTrack.tilesetInfo.coordSystem;
                            console.log(viewportTrack.tilesetInfo.coordSystem);
                        } catch(e) {

                        }
                        if(!assembly) return;

                        return (
                            <div className="col-tools-brush" key={i}
                                style={{
                                    // TODO: bar going out of the reange.
                                    left: startX,
                                    width: endX - startX
                                }}>
                                <div className={"chw-button-sm-container-horizontal"}
                                    style={{
                                        right: "4px",
                                        top: "2px",
                                        opacity: 1,
                                        background: "none",
                                        boxShadow: "none",
                                        color: "gray"
                                    }}>
                                    <svg className={`chw-button-sm-no-hover`}
                                        onClick={() => onRequestIntervalTFs({
                                            assembly,
                                            chrStartName,
                                            chrStartPos,
                                            chrEndName,
                                            chrEndPos
                                        })}
                                        viewBox={SEARCH.viewBox}>
                                        <title>Search bind TFs from Cistrome DB</title>
                                        <path d={SEARCH.path} fill="currentColor"/>
                                    </svg>
                                    <svg className={`chw-button-sm-no-hover`}
                                        onClick={() => onViewportRemove(viewportTrack.id)} 
                                        viewBox={CLOSE.viewBox}>
                                        <title>Remove viewport projection track</title>
                                        <path d={CLOSE.path} fill="currentColor"/>
                                    </svg>
                                </div>
                            </div>
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