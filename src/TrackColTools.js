import React from 'react';

import './TrackColTools.scss';

/**
 * Component for rendering genome interval selection tools.
 * @prop {number} viewY
 * @prop {number} viewHeight
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {(object|null)} combinedTrack The `combined` track, parent of the `horizontal-multivec` track.
 * @prop {object[]} siblingTracks An array of `viewport-projection-horizontal` track objects, which
 *                                are siblings of `multivecTrack` (children of the same `combined` track).
 * @prop {string} colToolsPosition The value of the `colToolsPosition` option.
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackColTools(props) {

    const {
        trackX, trackY, 
        trackWidth, trackHeight,
        trackAssembly,
        absGenomeScale,
        colToolsPosition,
        onSelectGenomicInterval
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

    const [g0, g1] = absGenomeScale.domain();

    function onButtonClick() {
        const gRange = (g1 - g0);
        const midpt = g0 + (gRange/2);
        onSelectGenomicInterval([ midpt - (gRange/10) , midpt + (gRange/10) ], undefined);
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
                <button 
                    onClick={onButtonClick}
                    style={{
                        marginTop: (isTop ? (2*height/3) : 2),
                    }}
                >Select interval</button>
            </div>
        </div>
    );
}