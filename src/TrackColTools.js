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
 */
export default function TrackColTools(props) {

    const {
        trackX, trackY, 
        trackWidth, trackHeight,
        trackAssembly,
        combinedTrack,
        siblingTracks,
        colToolsPosition,
        onSelectGenomicInterval
    } = props;

    if(!trackAssembly) {
        console.warn("trackAssembly is null");
        return null;
    }

    const left = trackX;
    const width = trackWidth;
    const height = 40;

    let top;
    if(colToolsPosition === "top") {
        top = trackY - height;
    } else if(colToolsPosition === "bottom") {
        top = trackY + trackHeight;
    }
    
    return (
        <div
            className="cistrome-hgw-child"
            style={{
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height}px`
            }}
        >
            <div className="col-tools">
                {!combinedTrack ? (
                    <button onClick={onSelectGenomicInterval}>Select genomic interval</button>
                ) : (
                    <div className="col-tools-selection-info">
                        {siblingTracks.map((siblingTrack, i) => (
                            <TrackColSelectionInfo
                                key={i}
                                trackAssembly={trackAssembly}
                                projectionTrack={siblingTrack}
                            />
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
};