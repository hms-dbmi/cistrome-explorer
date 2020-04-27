import React, { useRef, useEffect, useCallback } from 'react';



export default function TrackRowZoomOverlay(props) {
    const {
        trackX, trackY,
        trackWidth, trackHeight,
        onZoomRows,
        isWheelListening,
        numRowsTotal, numRowsSelected
    } = props;

    const top = trackY;
    const left = trackX;
    const height = trackHeight;
    const width = trackWidth;

    const overlayRef = useRef();

    useEffect(() => {
        const wheelHandler = (event) => {
            const { deltaMode, deltaY, layerY } = event;
            if(numRowsSelected !== undefined && numRowsTotal !== undefined && trackHeight !== undefined) {

                const zoomLevel = Math.max(0, Math.min(numRowsTotal, numRowsSelected - Math.floor(deltaY)));
                const zoomCenter = Math.floor((layerY / trackHeight) * numRowsSelected);
                console.log("zoomLevel", zoomLevel, "zoomCenter", zoomCenter);
                onZoomRows(zoomLevel, zoomCenter);
            }
        };
        if(isWheelListening) {
            overlayRef.current.addEventListener("wheel", wheelHandler);
        }
        return () => {
            overlayRef.current.removeEventListener("wheel", wheelHandler);
        }
    }, [overlayRef, isWheelListening, onZoomRows, trackHeight, numRowsTotal, numRowsSelected]);



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