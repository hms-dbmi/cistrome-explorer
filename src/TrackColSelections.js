import React from 'react';

import TrackColSelection from './TrackColSelection.js';


/**
 * Component for rendering multiple genome interval selections.
 * @prop {number} viewY
 * @prop {number} viewHeight
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackColSelections(props) {

    const {
        viewY,
        viewHeight,
        trackX, trackY,
        trackWidth, trackHeight,
        trackAssembly,
        absGenomeScale,
        colSelect,
        onSelectGenomicInterval,
        drawRegister
    } = props;

    if(!trackAssembly) {
        console.warn("trackAssembly is null");
        return null;
    }

    console.log("TrackColSelections.render");
    return (
        <div>
            {colSelect ? colSelect.map((interval, i) => (
                <TrackColSelection
                    key={i}
                    interval={interval}
                    viewY={viewY}
                    viewHeight={viewHeight}
                    trackX={trackX}
                    trackY={trackY}
                    trackWidth={trackWidth}
                    trackHeight={trackHeight}
                    trackAssembly={trackAssembly}
                    absGenomeScale={absGenomeScale}
                    onSelectGenomicInterval={(interval) => {
                        onSelectGenomicInterval(interval, i)
                    }}
                    drawRegister={drawRegister}
                />
            )) : null}
        </div>
    );
}