import React from "react";

/**
 * 
 * @param {*} props 
 */
export default function TrackResizer(props) {
    const {
        top, left, 
        isVisible
    } = props;

    return (
        <div
            className="visualization-resizer"
            onMouseDown={() => setOnResize(true)}
            style={{
                top: `${top + (height + 10) / 2.0}px`,
                left: `${isLeft ? 0 : width - 4}px`,
                height: `${10}px`,
                width: `${4}px`,
                visibility: mouseX !== null ? "visible" : "hidden"
            }}
        />
    );
} 