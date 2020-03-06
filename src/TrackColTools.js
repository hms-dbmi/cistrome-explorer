import React from 'react';

import TrackColSelectionInfo from './TrackColSelectionInfo.js';

import './TrackColTools.scss';

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
export default function TrackColTools(props) {

    const {
        trackX, trackY, 
        trackWidth, trackHeight,
        trackAssembly,
        combinedTrack,
        siblingTracks,
        colToolsPosition,
        onSelectGenomicInterval,
        onRequestIntervalTFs,
        drawRegister
    } = props;

    if(!trackAssembly) {
        console.warn("trackAssembly is null");
        return null;
    }

    const isTop = (colToolsPosition === "top");
    const left = trackX;
    const width = trackWidth;
    const height = 70;

    let top;
    if(colToolsPosition === "top") {
        top = trackY - height;
    } else if(colToolsPosition === "bottom") {
        top = trackY + trackHeight;
    }
    
    return (
        <div
            style={{
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height}px`,
                pointerEvents: 'none',
                position: "absolute"
            }}
        >
            <div className="col-tools">
                {!combinedTrack ? (
                    <button 
                        className="col-tools-target"
                        onClick={onSelectGenomicInterval}
                        style={{
                            marginTop: (isTop ? (2 * height / 3) : 2)
                        }}
                    >Select current interval</button>
                ) : (
                    <div className="col-tools-selection-info">
                        {siblingTracks.map((siblingTrack, i) => (
                            <TrackColSelectionInfo
                                key={i}
                                width={trackWidth}
                                height={height}
                                colToolsPosition={colToolsPosition}
                                trackAssembly={trackAssembly}
                                projectionTrack={siblingTrack}
                                onRequestIntervalTFs={onRequestIntervalTFs}
                                drawRegister={drawRegister}
                            />
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
};