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
        colSelect,
        trackAssembly,
        drawRegister
    } = props;

    if(!trackAssembly) {
        console.warn("trackAssembly is null");
        return null;
    }

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
                    drawRegister={drawRegister}
                />
            )) : null}
        </div>
    );
}