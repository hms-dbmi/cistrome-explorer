import React, { useEffect, useRef, useMemo, useCallback } from "react";
import d3 from "./utils/d3.js";

import TrackColSelectionInfo from "./TrackColSelectionInfo.js";

const HIGLASS_VIEWHEADER_HEIGHT = 24 + 4 + 4;

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
        absGenomeScale,
        onSelectGenomicInterval,
        drawRegister
    } = props;

    const gRef = useRef();
    const brush = useMemo(() => d3.brushX());

    const onBrushEnd = useCallback(() => {
        console.log("brushEnd");
        const event = d3.event;
        const viewportSelection = event.selection;
        const genomicSelection = viewportSelection.map(absGenomeScale.invert);

        onSelectGenomicInterval(genomicSelection);
    }, [absGenomeScale, onSelectGenomicInterval]);

    const [g0, g1] = absGenomeScale.domain();
    const [x0, x1] = absGenomeScale.range();


    useEffect(() => {
        brush
            .extent([
                [x0, viewY], 
                [x1, viewY + viewHeight]
            ]);
    }, [absGenomeScale, viewY, viewHeight]);


    useEffect(() => {
        const g = d3.select(gRef.current);
        
        g.call(brush);
        brush.on("end", null);
        
        // Check how much of the interval is within the current viewport.
        let viewportInterval;
        if(interval[0] >= g0 && interval[1] <= g1) {
            // Entire interval is within viewport.
            viewportInterval = interval.map(absGenomeScale);
        } else if(interval[0] < g0 && interval[1] > g0) {
            // Interval hangs off the left side of viewport.
            viewportInterval = [ x0, absGenomeScale(interval[1]) ];
        } else if(interval[0] < g1 && interval[1] > g1) {
            // Interval hangs off the right side of viewport.
            viewportInterval = [ absGenomeScale(interval[0]), x1 ];
        } else {
            // Interval is outside of current viewport.
            return;
        }

        g.call(brush.move, viewportInterval);
        brush.on("end", onBrushEnd);

        // Turn off the ability to select new regions for this brush
        // See https://github.com/higlass/higlass/blob/develop/app/scripts/ViewportTrackerHorizontal.js#L34
        g.selectAll('.overlay')
            .style('pointer-events', 'none');

    }, [interval, absGenomeScale]);

    console.log("TrackColSelection.render");
    return (
        <div>
            <svg
                style={{
                    position: 'absolute',
                    top: `${HIGLASS_VIEWHEADER_HEIGHT + viewY}px`,
                    left: `${trackX}px`,
                    height: `${viewHeight}px`,
                    width: `${trackWidth}px`,
                    pointerEvents: 'none',
                }}
            >
                <g 
                    className="brush" 
                    ref={gRef} 
                />
            </svg>
        </div>
    );
}