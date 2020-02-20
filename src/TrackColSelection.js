import React, { useEffect, useRef } from "react";
import d3 from "./utils/d3.js";

import TrackColSelectionInfo from "./TrackColSelectionInfo.js";

/**
 * Component for rendering genome interval selection brushing elements.
 * @prop {number} viewY
 * @prop {number} viewHeight
 */
export default function TrackColSelection(props) {
    const {
        interval,
        viewY,
        viewHeight,
        trackX, trackY, 
        trackWidth, trackHeight,
        trackAssembly,
        drawRegister
    } = props;

    const hgViewHeaderHeight = 24 + 4 + 4;
    const [intervalStart, intervalEnd] = interval;
    const intervalWidth = intervalEnd - intervalStart;

    const svgRef = useRef();

    useEffect(() => {
        const svg = svgRef.current;
        svg.style.left = `${intervalStart}px`;
        svg.style.width = `${intervalWidth}px`;
        
    });

    return (
        <div>
            <svg 
                ref={svgRef}
                style={{
                    position: 'absolute',
                    top: `${hgViewHeaderHeight + viewY}px`,
                    height: `${viewHeight}px`,
                    backgroundColor: 'blue',
                    opacity: 0.5,
                    pointerEvents: 'none',
                }}
            />
        </div>
    );
}