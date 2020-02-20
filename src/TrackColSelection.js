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

    function brushed(e) {
        console.log(e);
    }

    useEffect(() => {
        const svg = d3.select(svgRef.current);
 
        const g = svg.append("g")
            .attr("class", "brush");

        const brush = d3.brushX()
            .extent([[trackX, viewY], [trackX + trackWidth, viewY + viewHeight]]);
        
        g.call(brush);
        brush.on("brush", null);
        g.call(brush.move, interval);
        brush.on("brush", brushed);

        // turn off the ability to select new regions for this brush
        g.selectAll('.overlay')
            .style('pointer-events', 'none');
            
    }, [svgRef]);

    return (
        <div>
            <svg
                ref={svgRef}
                style={{
                    position: 'absolute',
                    top: `${hgViewHeaderHeight + viewY}px`,
                    left: `${trackX}px`,
                    height: `${viewHeight}px`,
                    width: `${trackWidth}px`,
                    pointerEvents: 'none',
                }}
            />
        </div>
    );
}