import React, { useRef, useEffect, useCallback } from 'react';
import debounce from 'lodash/debounce';

export default function TrackRowZoomOverlay(props) {
    const {
        trackX, trackY,
        trackWidth, trackHeight,
        onZoomRows,
        isWheelListening
    } = props;

    const top = trackY;
    const left = trackX;
    const height = trackHeight;
    const width = trackWidth;

    const overlayRef = useRef();

    useEffect(() => {
        const wheelHandler = debounce((event) => {
            const { deltaMode, deltaY, layerY } = event;
            onZoomRows(layerY / trackHeight, deltaY, deltaMode);
        }, 50);
        if(isWheelListening) {
            overlayRef.current.addEventListener("wheel", wheelHandler);
        }
        return () => {
            overlayRef.current.removeEventListener("wheel", wheelHandler);
        }
    }, [overlayRef, isWheelListening, onZoomRows, trackHeight]);



    const style = (isWheelListening ? {
        position: 'absolute',
        top: `${top}px`,
        left: `${left}px`,
        width: `${width}px`,
        height: `${height}px`,
        backgroundColor: 'rgba(138,43,226, 0.4)',
        display: 'block',
        lineHeight: 1
    } : {
        display: 'none'
    });

    return (<div
        ref={overlayRef}
        style={style}
    />);
}