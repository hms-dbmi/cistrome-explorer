import React, { useRef, useState, useEffect, useCallback } from 'react';
import d3 from './utils/d3.js';
import TrackColSelectionInfo from './TrackColSelectionInfo.js';

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

    // TODO: All hooks must be above this return statement, since they need to be executed in the same order.
    if(!viewBoundingBox) {
        // The view info has not yet loaded.
        return null;
    }

    // if(!trackAssembly) {
    //     console.warn("trackAssembly is null");
    //     return null;
    // }

    const divRef = useRef();
    const [mouseX, setMouseX] = useState(null);
    const [mouseChrX, setMouseChrX] = useState(null);
    const [brushStartX, setBrushStartX] = useState(null);
    const [brushWidth, setBrushWidth] = useState(0);
    const dragX = useRef(null);

    // TODO: Remove this.
    const isTop = (colToolsPosition === "top");
    
    const { top, left, width, height } = viewBoundingBox;

    const brushBarTop = height + 4;
    const brushBarHeight = 30;
    
    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const div = divRef.current;
        const event = d3.event;
        dragX.current = event.sourceEvent.clientX;

        const [mouseX, mouseY] = d3.mouse(div);
        setBrushStartX(mouseX);
        setBrushWidth(0);
    }, [dragX])

    const dragged = useCallback(() => {
        const event = d3.event;
        const diff = event.sourceEvent.clientX - dragX.current;
        let newWidth = brushWidth + diff;
        setBrushWidth(newWidth);
    }, [dragX]);

    const ended = useCallback(() => {
        dragX.current = null;
    }, [dragX])

    // Detect drag events for the resize element.
    useEffect(() => {
        const div = divRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(div).call(drag);

        return () => d3.select(resizer).on(".drag", null);
    }, [divRef, started, dragged, ended]);

    useEffect(() => {
        const div = divRef.current;

        d3.select(div).on("mousemove", () => {
            const [mouseX, mouseY] = d3.mouse(div);

            // TODO: Remove this.
            // console.log(multivecTrack);

            // TODO: Change this to the actual chr position.
            // const xAbsPos = multivecTrack._xScale.invert(mouseX);

            setMouseX(mouseX);
            // setMouseChrX(xAbsPos);
        });
        // TODO: uncomment when ready to release
        // d3.select(div).on("mouseleave", () => {
        //     setMouseX(null);
        //     setMouseChrX(null);
        // });
    }, [divRef]);

    return (
        <div style={{
            position: "absolute",
            top: `${top}px`,
            left: `${left}px`,
            width: `${width}px`, 
            height: `${height}px`
        }}>
            <div className="col-tools" ref={divRef}
                style={{
                    top: `${brushBarTop}px`,
                    height: `${brushBarHeight}px`,
                }}>
                {mouseChrX ? 
                    <div className="col-tools-chr-info" style={{
                        left: `${mouseX + 16}px`,
                        lineHeight: `${brushBarHeight}px`
                    }}>
                        {(+mouseChrX.toFixed()).toLocaleString("en")}
                    </div>
                    : null}
                {/* TODO: Remove Code Below when ready to release */}
                <button 
                    className="col-tools-target"
                    onClick={onSelectGenomicInterval}
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
                ) : null}
                {/* TODO: remove above. */}
            </div>
            <div className="col-tools-hg-wrapper">
                {mouseX ? 
                    <div className="col-tools-column-info" style={{
                        left: `${mouseX}px`
                    }}/>
                    : null}
                {brushWidth !== 0 ?
                    <div className="col-tools-brush" style={{
                        left: brushStartX,
                        width: brushWidth
                    }}/>
                    : null}
            </div>
        </div>
    );
};