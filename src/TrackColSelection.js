import React, { useEffect, useRef } from "react";
import d3 from "./utils/d3.js";

/**
 * Component for rendering genome interval selection brushing elements.
 * @prop {number} viewY
 * @prop {number} viewHeight
 */
export default function TrackColSelection(props) {
    const {
        viewY,
        viewHeight,
        interval
    } = props;

    const hgViewHeaderHeight = 24 + 4 + 4;
    const [intervalStart, intervalEnd] = interval;
    const intervalWidth = intervalEnd - intervalStart;

    return (
        <div
            style={{
                position: 'absolute',
                top: `${hgViewHeaderHeight + viewY}px`,
                height: `${viewHeight}px`,
                left: `${intervalStart}px`,
                width: `${intervalWidth}px`,
                backgroundColor: 'blue',
                opacity: 0.5,
                pointerEvents: 'none',
            }}
        />
    );
}